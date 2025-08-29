#region modules
from lark import Transformer, Lark, Token, v_args
from collections import OrderedDict
from typing import Dict, List, Any
#endregion

#region variables
#endregion

#region functions
def _fmt(x: Any) -> str:
    if isinstance(x, (int, float)):
        return str(int(x)) if float(x).is_integer() else f"{x:.10g}"
    return str(x)
#endregion

#region classes
class AbacusInputTransform(Transformer):
    def __init__(self):
        super().__init__(visit_tokens=True)
        self.dict = {}

    def NAME(self, args): return args.value
    def VALUE(self, args): return args.value.strip(' ') 

    def pair(self, args):
        return [args[0], args[1]]
    
    def start(self, args):
        args = list(filter(lambda item: isinstance(item, list), args))
        for pair in args:
            key = pair[0]
            value = pair[1]

            self.dict.update({key: value})

        return self.dict

@v_args(inline=True)
class AbacusStruTransform(Transformer):
    # ---- tokens ----
    def NUMBER(self, t):  # handles 1, 1.0, 1e-3, 1D-3, etc.
        return float(t.value.replace('D', 'e').replace('d', 'e'))
    def SYMBOL(self, t): return t.value
    def NAME(self, t):   return t.value
    def NEWLINE(self, _): return None  # drop

    # ---- small helpers ----
    def _nz(self, xs):  # non-None
        return [x for x in xs if x is not None]

    # ---- line rules ----
    def species_line(self, sym, mass, name, *_):
        return [sym, mass, name]

    def orbital_line(self, name, *_):
        return [name]

    def vectors_line(self, *nums):
        nums = self._nz(nums)
        return nums  # [x, y, z]

    def position_vector(self, *nums):
        nums = self._nz(nums)
        if len(nums) != 6:
            raise ValueError(f"position_vector expects 6 numbers, got {len(nums)}: {nums}")
        return nums  # [x, y, z, fx, fy, fz]

    def position_entry(self, *kids):
        kids = self._nz(kids)
        sym, mm, count, *vecs = kids
        count = int(round(count))
        if len(vecs) < count:
            raise ValueError(f"{sym}: need {count} vectors, got {len(vecs)}")
        rows = []
        for v in vecs[:count]:
            x, y, z, fx, fy, fz = v
            rows.append([sym, x, y, z, fx, fy, fz, mm])
        return rows

    # ---- block rules -> (KEY, value) ----
    def atomic_species(self, *_kids):
        lines = [c for c in _kids if isinstance(c, list)]
        return ("ATOMIC_SPECIES", lines)

    def numerical_orbital(self, *_kids):
        lines = [c for c in _kids if isinstance(c, list)]
        return ("NUMERICAL_ORBITAL", lines)

    def lattice_constant(self, *_kids):
        nums = [c for c in _kids if isinstance(c, (int, float))]
        if len(nums) != 1:
            raise ValueError(f"LATTICE_CONSTANT expects 1 number, got {nums}")
        return ("LATTICE_CONSTANT", [[nums[0]]])

    def lattice_vectors(self, *_kids):
        lines = [c for c in _kids if isinstance(c, list)]
        return ("LATTICE_VECTORS", lines)

    def atomic_positions(self, *_kids):
        kids = self._nz(_kids)
        # find coordinate type (e.g., 'Direct'); grammar has a literal "Direct"
        coord_type = None
        for k in kids:
            val = k.value if isinstance(k, Token) else k
            if isinstance(val, str) and val.strip():
                coord_type = val.strip()
                break
        data = []
        for k in kids:
            if isinstance(k, list):
                # position_entry returns a list of rows
                if k and isinstance(k[0], list):
                    data.extend(k)
        return ("ATOMIC_POSITIONS", {"type": coord_type or "Direct", "data": data})

    # ---- aggregation ----
    def item(self, kv):
        return kv  # (key, value)

    def items(self, *kvs):
        d = {}
        for k, v in kvs:
            d[k] = v
        return d

    def start(self, *children):
        for c in children:
            if isinstance(c, dict):
                return c
        return {}

class AbacusKptTransform(Transformer):
    def __init__(self):
        super().__init__(visit_tokens=True)
        self.dict = {}

    def NUMBER(self, args): return args.value
    def TYPE(self, args): return args.value

    def gamma_line(self, args):
        return args[:-1]
    
    def kpt_line(self, args):
        return args[:-1]
    
    def start(self, args):
        nkpt = args[1]
        kpt_type = args[3]
        data_list = args[5:]

        self.dict.update({'nkpt': nkpt})
        self.dict.update({'type': kpt_type})
        self.dict.update({'data': data_list})

        return self.dict

class AbacusInputGrammar:
    grammar = r'''
%import common.NEWLINE
%import common.WS_INLINE
%ignore WS_INLINE

NAME: /[A-Za-z_][A-Za-z0-9_\-\.]*/
VALUE: /.+/

start: "INPUT_PARAMETERS" NEWLINE pair+ 

pair: NAME VALUE NEWLINE
'''
    transform = AbacusInputTransform()

    def __init__(self):
        self.parser = Lark(self.grammar, parser='lalr')

    def read(self, text: str) -> dict:
        tree = self.parser.parse(text)
        return self.transform.transform(tree)

    def write(self, data: dict) -> str:
        output = f'INPUT_PARAMETERS\n'

        for key, value in data.items():
            output += f'{key} {value}\n'

        return output 

class AbacusStruGrammar:
    grammar = r'''
%import common.NEWLINE
%import common.WS_INLINE
%ignore WS_INLINE

NUMBER: /[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?/
SYMBOL: /[A-Z][a-z]?/
NAME: /[A-Za-z_][A-Za-z0-9_\-\.]*/

start: NEWLINE* items NEWLINE*

items: item+

item: atomic_species | numerical_orbital | lattice_constant | lattice_vectors | atomic_positions

atomic_species: "ATOMIC_SPECIES" NEWLINE species_line+
species_line: SYMBOL NUMBER NAME NEWLINE

numerical_orbital: "NUMERICAL_ORBITAL" NEWLINE orbital_line+ 
orbital_line: NAME NEWLINE

lattice_constant: "LATTICE_CONSTANT" NEWLINE NUMBER NEWLINE

lattice_vectors: "LATTICE_VECTORS" NEWLINE vectors_line+
vectors_line: NUMBER NUMBER NUMBER NEWLINE

atomic_positions: "ATOMIC_POSITIONS" NEWLINE ("Direct" | "Cartesian") NEWLINE position_entry+ 
position_entry: SYMBOL NEWLINE NUMBER NEWLINE NUMBER NEWLINE position_vector+
position_vector: NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NEWLINE
'''
    transform = AbacusStruTransform()

    def __init__(self):
        self.parser = Lark(self.grammar, parser='lalr')

    def read(self, text: str) -> dict:
        tree = self.parser.parse(text)
        return self.transform.transform(tree)

    def write(self, data: Dict) -> str:
        parts: List[str] = []

        # --- ATOMIC_SPECIES ---
        if "ATOMIC_SPECIES" in data:
            lines = ["ATOMIC_SPECIES"]
            for sym, mass, pseudo in data["ATOMIC_SPECIES"]:
                lines.append(f"{sym} {_fmt(mass)}  {pseudo}")
            parts.append("\n".join(lines))

        # --- NUMERICAL_ORBITAL ---
        if "NUMERICAL_ORBITAL" in data:
            lines = ["", "NUMERICAL_ORBITAL"]
            for (name,) in data["NUMERICAL_ORBITAL"]:
                lines.append(name)
            parts.append("\n".join(lines))

        # --- LATTICE_CONSTANT ---
        if "LATTICE_CONSTANT" in data:
            val = data["LATTICE_CONSTANT"][0][0]
            parts.append("\n".join(["", "LATTICE_CONSTANT", _fmt(val)]))

        # --- LATTICE_VECTORS ---
        if "LATTICE_VECTORS" in data:
            lines = ["", "LATTICE_VECTORS"]
            for x, y, z in data["LATTICE_VECTORS"]:
                lines.append(f"{_fmt(x)} {_fmt(y)} {_fmt(z)}")
            parts.append("\n".join(lines))

        # --- ATOMIC_POSITIONS ---
        if "ATOMIC_POSITIONS" in data:
            ap = data["ATOMIC_POSITIONS"]
            ctype = ap.get("type", "Direct")
            rows = ap.get("data", [])

            # Group by (symbol, magnetic_moment) in first-seen order
            groups: "OrderedDict[tuple, List[List[float]]]" = OrderedDict()
            for sym, x, y, z, fx, fy, fz, mm in rows:
                key = (sym, mm)
                groups.setdefault(key, []).append([x, y, z, fx, fy, fz])

            lines = ["", "ATOMIC_POSITIONS", ctype]
            for (sym, mm), vecs in groups.items():
                lines.append(sym)
                lines.append(_fmt(mm))
                lines.append(str(len(vecs)))
                for x, y, z, fx, fy, fz in vecs:
                    lines.append(
                        f"{_fmt(x)}  {_fmt(y)}  {_fmt(z)}  {_fmt(fx)} {_fmt(fy)} {_fmt(fz)}"
                    )
            parts.append("\n".join(lines))

        # Join all blocks with a trailing newline like typical INPUTs
        out = "\n".join(parts).lstrip("\n")
        return out + ("\n" if not out.endswith("\n") else "")

class AbacusKptGrammar:
    grammar = r'''
%import common.NEWLINE
%import common.WS_INLINE
%ignore WS_INLINE

NUMBER: /[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?/
TYPE: "Line" | "Direct" | "Gamma"

start: "K_POINTS" NEWLINE NUMBER NEWLINE TYPE NEWLINE (gamma_line | kpt_line)+ 

gamma_line: NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NEWLINE

kpt_line: NUMBER NUMBER NUMBER NUMBER NEWLINE
'''
    transform = AbacusKptTransform()

    def __init__(self):
        self.parser = Lark(self.grammar, parser='lalr')

    def read(self, text: str) -> dict:
        tree = self.parser.parse(text)
        return self.transform.transform(tree)

    def write(self, data: dict) -> str:
        output = f'K_POINTS\n'
        output += f'{data['nkpt']}\n'
        output += f'{data['type']}\n'

        for list_item in data['data']:
            list_item_str = ' '.join(list_item)
            output += f'{list_item_str}\n'

        return output 

#endregion