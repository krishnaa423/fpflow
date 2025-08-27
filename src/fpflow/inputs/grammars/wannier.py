# region modules
from lark import Lark, Transformer, Token
import math
# endregion

#region: variables
#endregion

#region: functions
#endregion

# region classes
class WannierTransform(Transformer):
    """Transforms the parse tree into the requested dict schema."""

    # ----- token converters -----
    def NUMBER(self, t: Token):
        s = str(t).replace('D', 'E').replace('d', 'e')
        v = float(s)
        if math.isfinite(v) and abs(v - round(v)) < 1e-12:
            return int(round(v))
        return v

    def LOGICAL(self, t: Token):
        return str(t).lower() == ".true."

    def VALUE_STR(self, t: Token):
        return str(t)[1:-1]  # strip single quotes

    def NAME(self, t: Token):
        return str(t)

    def SYMBOL(self, t: Token):
        return str(t)

    def SPECIAL_POINT(self, t: Token):
        return str(t)

    def UNIT(self, t: Token):
        return str(t)

    # ----- small lines -> lists -----
    def cell_line(self, items):
        return [items[0], items[1], items[2]]

    def atoms_line(self, items):
        return [items[0], items[1], items[2], items[3]]

    def kpoints_line(self, items):
        return [items[0], items[1], items[2]]

    def kpath_line(self, items):
        return [items[0], items[1], items[2], items[3],
                items[4], items[5], items[6], items[7]]

    # ----- helpers for values -----
    def number_list(self, items):
        return list(items)

    def value(self, items):
        return items[0]

    def unit(self, items):
        return items[0]  # already converted by UNIT()

    # ----- blocks -> dict entries -----
    def unit_cell_cart(self, items):
        unit = items[1]
        data = [row for row in items[3:]]
        return {"unit_cell_cart": {"unit": unit, "data": data}}

    def atoms_frac(self, items):
        return {"atoms_frac": {"unit": "frac", "data": list(items[1:])}}

    def kpoints(self, items):
        return {"kpoints": list(items[1:])}

    def kpath(self, items):
        return {"kpoint_path": list(items[1:])}

    # ----- pair -> dict entry -----
    def pair(self, items):
        key, val = items
        return {key: val}

    # ----- wrapper rules -----
    def block(self, items):
        # Unwrap the single dict produced by a block rule
        return items[0]

    def statement(self, items):
        # Should already be a dict
        return items[0]

    # ----- start: defensively merge dicts -----
    def start(self, items):
        out = {}
        for piece in items:
            if isinstance(piece, dict):
                out.update(piece)
                continue
            # Some parser/transformer inlining could wrap a dict in a list
            if isinstance(piece, (list, tuple)):
                for it in piece:
                    if isinstance(it, dict):
                        out.update(it)
                continue
            # Last resort: ignore anything that isn't a dict
        return out

class WannierGrammar:
    grammar: str = r"""
%import common.NEWLINE
%import common.WS_INLINE
%ignore WS_INLINE

COMMENT: /![^\n]*/
%ignore COMMENT

NAME: /[A-Za-z_][A-Za-z0-9_\.]*/
LOGICAL: /\.(?:true|false)\./i
VALUE_STR: "'" /[^']*/ "'"
SYMBOL: /[A-Z][a-z]?/
SPECIAL_POINT: /[A-Z]/
NUMBER: /[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?/

start: NEWLINE* statement (NEWLINE+ statement)* NEWLINE*

statement: pair | block

pair: NAME "=" value
value: VALUE_STR | NUMBER | number_list | NAME | LOGICAL
number_list: NUMBER (NUMBER)+

block: unit_cell_cart | atoms_frac | kpoints | kpath

UNIT: "bohr" | "angstrom"
unit: UNIT

unit_cell_cart: "begin" "unit_cell_cart" NEWLINE unit NEWLINE cell_line+ "end" "unit_cell_cart"
cell_line: NUMBER NUMBER NUMBER NEWLINE

atoms_frac: "begin" "atoms_frac" NEWLINE atoms_line+ "end" "atoms_frac"
atoms_line: SYMBOL NUMBER NUMBER NUMBER NEWLINE

kpoints: "begin" "kpoints" NEWLINE kpoints_line+ "end" "kpoints"
kpoints_line: NUMBER NUMBER NUMBER NEWLINE

kpath: "begin" "kpoint_path" NEWLINE kpath_line+ "end" "kpoint_path"
kpath_line: SPECIAL_POINT NUMBER NUMBER NUMBER SPECIAL_POINT NUMBER NUMBER NUMBER NEWLINE
"""

    transform: WannierTransform = WannierTransform()

    def __init__(self):
        self.parser = Lark(self.grammar, parser='lalr')

    def read(self, text: str) -> dict:
        tree = self.parser.parse(text)
        return self.transform.transform(tree)

    def write(self, data: dict) -> str:
        """Serialize dict -> Wannier90 input text."""
        lines = []

        def fmt_scalar(v):
            if isinstance(v, bool):
                return ".true." if v else ".false."
            if isinstance(v, (int, float)):
                return str(v)
            return str(v)

        def fmt_row(row):
            return " ".join(fmt_scalar(x) for x in row)

        def fmt_kv(key, val):
            if isinstance(val, (list, tuple)):
                return f"{key} = " + " ".join(fmt_scalar(x) for x in val)
            else:
                return f"{key} = {fmt_scalar(val)}"

        block_keys = {"unit_cell_cart", "atoms_frac", "kpoints", "kpoint_path"}

        # 1) key-value pairs (everything not in block_keys)
        for key, val in data.items():
            if key in block_keys:
                continue
            lines.append(fmt_kv(key, val))

        # 2) blocks

        # unit_cell_cart
        if "unit_cell_cart" in data:
            blk = data["unit_cell_cart"] or {}
            unit = blk.get("unit", "bohr")
            cell = blk.get("data", [])
            lines.append("\nbegin unit_cell_cart")
            lines.append(unit)
            for row in cell:
                lines.append(fmt_row(row))
            lines.append("end unit_cell_cart\n")

        # atoms_frac
        if "atoms_frac" in data:
            blk = data["atoms_frac"] or {}
            rows = blk.get("data", [])
            lines.append("\nbegin atoms_frac")
            for row in rows:
                if not row:
                    continue
                sym = row[0]
                coords = row[1:]
                lines.append(" ".join([str(sym)] + [fmt_scalar(x) for x in coords]))
            lines.append("end atoms_frac\n")

        # kpoints
        if "kpoints" in data:
            rows = data["kpoints"] or []
            lines.append("\nbegin kpoints")
            for row in rows:
                lines.append(fmt_row(row))
            lines.append("end kpoints\n")

        # kpoint_path
        if "kpoint_path" in data:
            rows = data["kpoint_path"] or []
            lines.append("\nbegin kpoint_path")
            for row in rows:
                if len(row) == 8:
                    P, kx, ky, kz, Q, qx, qy, qz = row
                    lines.append(
                        f"{P} {fmt_scalar(kx)} {fmt_scalar(ky)} {fmt_scalar(kz)} "
                        f"{Q} {fmt_scalar(qx)} {fmt_scalar(qy)} {fmt_scalar(qz)}"
                    )
                else:
                    # fallback: write whatever is present
                    lines.append(" ".join(str(x) for x in row))
            lines.append("end kpoint_path\n")

        return "\n".join(lines) + "\n"

# endregion
