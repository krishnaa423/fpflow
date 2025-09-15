#region modules
from __future__ import annotations
from lark import Lark, Transformer, Token
import re 
from typing import Any, Dict, List
from fpflow.io.logging import get_logger
#endregion

#region variables
logger = get_logger()
#endregion

#region functions
#endregion

#region classes
class QeTransform(Transformer):
    """
    Transforms a QE input parsed with the user's grammar into:

    {
      "control": {...}, "system": {...}, "electrons": {...}, "ions": {...}, "cell": {...},
      "atomic_species": [ [sym, mass, file], ... ],
      "atomic_positions": [ [sym, x, y, z], ... ],
      "cell_parameters": [ [a11,a12,a13], [a21,a22,a23], [a31,a32,a33] ],
      "k_points": {"type":"crystal","nkpt":"27","data":[[kx,ky,kz,w], ...]}
        # or {"type":"automatic","data":["nx","ny","nz","sx","sy","sz"]}
    }

    All final values are strings.
    """

    def __init__(self):
        super().__init__(visit_tokens=True)

    # --- token mappers ---
    def NAME(self, t): return t.value
    def FILENAME(self, t): return t.value
    def SYMBOL(self, t): return t.value
    def NUMBER(self, t): return t.value
    def INT(self, t): return t.value
    def VALUE_STR(self, t): return t.value[1:-1]
    def LOGICAL(self, t): return t.value.strip(".").lower()
    def NEWLINE(self, _): return None

    @staticmethod
    def _str(x): return x if isinstance(x, str) else str(x)
    
    def item(self, items): return items[0] if items else None

    def start(self, items):
        out = {}
        for it in items:
            if isinstance(it, dict):
                out.update(it)
        return out

    # --- namelist same as before ---
    def namelist(self, items):
        name, body = None, {}
        for x in items:
            if isinstance(x, str) and name is None:
                name = x.lower()
            elif isinstance(x, dict):
                body = x
        return {name: body}

    def nl_body(self, items):
        d = {}
        for kv in items:
            if isinstance(kv, tuple) and len(kv) == 2:
                k, v = kv
                d[self._str(k)] = self._str(v)
        return d

    def nl_entry(self, items):
        parts = [x for x in items if x is not None and not (isinstance(x, Token) and x.value == "=")]
        return (self._str(parts[0]), self._str(parts[1]))

    def nl_key(self, items):
        name = self._str(items[0])
        if len(items) > 1:
            idx = next((self._str(x) for x in items[1:] if isinstance(x, str)), None)
            if idx: return f"{name}({idx})"
        return name

    def nl_value(self, items): return self._str(items[0])

    # --- ATOMIC_SPECIES ---
    def atomic_species(self, items):
        rows = []
        for row in items:
            if isinstance(row, list):
                rows.append([self._str(x) for x in row])
        return {"atomic_species": rows}

    def species_line(self, items):
        return [self._str(items[0]), self._str(items[1]), self._str(items[2])]

    # --- ATOMIC_POSITIONS ---
    def unit(self, items): return self._str(items[-1])

    def atomic_positions(self, items):
        unit, rows = None, []
        for x in items:
            if isinstance(x, str) and unit is None:
                unit = x
            elif isinstance(x, list):
                rows.append([self._str(v) for v in x])
        return {"atomic_positions": {"unit": unit, "data": rows}}

    def position_line(self, items):
        return [self._str(x) for x in items if x is not None]

    # --- CELL_PARAMETERS ---
    def cell_parameters(self, items):
        unit, rows = None, []
        for x in items:
            if isinstance(x, str) and unit is None:
                unit = x
            elif isinstance(x, list):
                rows.append([self._str(v) for v in x])
        return {"cell_parameters": {"unit": unit, "data": rows}}

    def cell_row(self, items):
        return [self._str(x) for x in items if x is not None]

    # --- k-points same as before ---
    def k_points(self, items):
        kunit, payload = None, None
        for x in items:
            if isinstance(x, str) and kunit is None:
                kunit = x.lower()
            elif payload is None:
                payload = x
        if kunit is None:
            kunit = "crystal" if isinstance(payload, tuple) else "automatic"
        if kunit == "automatic":
            return {"k_points": {"type": "automatic", "data": [self._str(v) for v in payload]}}
        nkpt, pts = payload
        return {"k_points": {"type": "crystal",
                             "nkpt": self._str(nkpt),
                             "data": [[self._str(v) for v in row] for row in pts]}}

    def k_auto(self, items): return [self._str(v) for v in items if v is not None][:6]
    
    def k_crystal(self, items):
        nkpt = self._str(items[0])
        pts = [row for row in items[1:] if isinstance(row, list)]
        return (nkpt, pts)
    
    def kpt_line(self, items): return [self._str(x) for x in items if x is not None][:4]

class QeGrammar:
    grammar = r'''
%import common.WS_INLINE
%import common.NEWLINE
%import common.INT
%ignore WS_INLINE

COMMENT: /![^\n]*/
%ignore COMMENT

NAME: /[A-Za-z_][A-Za-z0-9_\.]*/
FILENAME: /[^\s]+/
LOGICAL: /\.(?:true|false)\./i
VALUE_STR: "'" /[^']*/ "'"
SYMBOL: /[A-Z][a-z]?/
NUMBER: /[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?/

start: item+ 

item: namelist | atomic_species | cell_parameters | atomic_positions | k_points 

namelist: "&" NAME NEWLINE nl_body "/" NEWLINE 
nl_body: nl_entry* 
nl_entry: nl_key "=" nl_value NEWLINE 
nl_key: NAME ("(" INT ")")?
nl_value: VALUE_STR | NUMBER | NAME | LOGICAL

atomic_species: "ATOMIC_SPECIES" NEWLINE species_line+
species_line: SYMBOL NUMBER FILENAME NEWLINE

unit: NAME | "(" NAME ")"

atomic_positions: "ATOMIC_POSITIONS" unit? NEWLINE position_line+
position_line: SYMBOL NUMBER NUMBER NUMBER NEWLINE

cell_parameters: "CELL_PARAMETERS" unit? NEWLINE cell_row cell_row cell_row
cell_row: NUMBER NUMBER NUMBER NEWLINE

k_points: "K_POINTS" unit? NEWLINE (k_auto | k_crystal)
k_auto: INT INT INT INT INT INT NEWLINE
k_crystal: INT NEWLINE kpt_line+
kpt_line: NUMBER NUMBER NUMBER NUMBER NEWLINE
'''
    transform: QeTransform = QeTransform()

    def __init__(self):
        self.parser = Lark(self.grammar, parser='lalr')

    def read(self, text: str) -> dict:
        tree = self.parser.parse(text)
        return self.transform.transform(tree)
    
    def write(self, data: dict) -> str:
        def is_number(s: str): return bool(re.fullmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?", s))
        def fmt_value(v: str) -> str:
            lv = v.lower()
            if lv == "true": return ".true."
            if lv == "false": return ".false."
            if is_number(v): return v
            return f"'{v.replace('\'','\'\'')}'"

        def write_namelist(name, body):
            if body is None: return ""
            lines = [f"&{name.upper()}"]
            for k, v in body.items():
                lines.append(f"{k} = {fmt_value(str(v))}")
            lines.append("/\n")
            return "\n".join(lines) + "\n"

        out = []

        # Namelists
        for nl in ("control","system","electrons","ions","cell"):
            if nl in data:
                body = data[nl]
                if isinstance(body, dict):
                    out.append(write_namelist(nl, body))
                elif body == {}:
                    out.append(f"&{nl.upper()}\n/\n")

        out.append("\n")

        # ATOMIC_SPECIES
        species = data.get("atomic_species")
        if species:
            lines = ["ATOMIC_SPECIES"]
            for sym, mass, fname in species:
                lines.append(f"{sym} {mass} {fname}")
            out.append("\n".join(lines) + "\n\n")

        # ATOMIC_POSITIONS
        apos = data.get("atomic_positions")
        if apos and isinstance(apos, dict):
            unit = apos.get("unit")
            rows = apos.get("data", [])
            header = "ATOMIC_POSITIONS" + (f" {unit}" if unit else "")
            lines = [header] + [" ".join(list(map(str, row))) for row in rows]
            out.append("\n".join(lines) + "\n\n")

        # CELL_PARAMETERS
        cellp = data.get("cell_parameters")
        if cellp and isinstance(cellp, dict):
            unit = cellp.get("unit")
            rows = cellp.get("data", [])
            header = "CELL_PARAMETERS" + (f" {unit}" if unit else "")
            lines = [header] + [" ".join(list(map(str, row))) for row in rows]
            out.append("\n".join(lines) + "\n\n")

        # K_POINTS
        kpts = data.get("k_points")
        if kpts:
            if kpts.get("type") == "automatic":
                arr = kpts.get("data", [])
                out.append("K_POINTS automatic\n" + " ".join(list(map(str, arr))) + "\n")
            elif kpts.get("type") == "crystal":
                nkpt = kpts.get("nkpt")
                rows = kpts.get("data", [])
                lines = ["K_POINTS crystal", str(nkpt)] + [" ".join(list(map(str, row))) for row in rows]
                out.append("\n".join(lines) + "\n")

        # OCCUPATIONS
        occupations = data.get("occupations", None)
        if occupations is not None:
            out.append("\nOCCUPATIONS\n")
            for occ in occupations:
                out.append(f"{occ}\n")

        text = "".join(out)
        if not text.endswith("\n"):
            text += "\n"
        return text

#endregion