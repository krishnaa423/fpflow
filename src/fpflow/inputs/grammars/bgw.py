#region modules
from lark import Lark, Transformer
import re
#endregion

#region functions
def _to_number(s: str):
    if re.match(r'^[+-]?\d+$', s):
        return int(s)
    return float(s)

def _fmt_number(x):
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        return f"{x:.10g}"
    return str(x)
#endregion

#region classes
class BgwTransform(Transformer):
    def file(self, items):
        out = {}
        for k, v in items:
            out[k] = v
        return out

    def pair(self, items):
        key = str(items[0])
        val = items[1]
        return (key, val)

    def flag(self, items):
        key = str(items[0])
        return (key, None)

    def block(self, items):
        name = str(items[0])
        rows = [row for row in items[1:]]
        return (name, rows)

    def row(self, items):
        return [_to_number(str(tok)) for tok in items]

    def number(self, items):
        return _to_number(str(items[0]))


class BgwGrammar:
    grammar: str = r"""
start: _NL* statement (_NL+ statement)* _NL*   -> file

?statement: pair
          | flag
          | block

pair: NAME value                     -> pair
flag: NAME                           -> flag

block: "begin"i NAME _NL row+ "end"i -> block

row: SIGNED_NUMBER+ _NL              -> row

value: SIGNED_NUMBER                 -> number

_NL: (NL)+

NAME: /[A-Za-z_][A-Za-z0-9_\.]*/

%import common.SIGNED_NUMBER
%import common.NEWLINE -> NL
%import common.WS_INLINE
%ignore WS_INLINE
"""

    transform: BgwTransform = BgwTransform()

    def __init__(self):
        self.parser = Lark(self.grammar, parser='lalr')

    def read(self, text: str) -> dict:
        tree = self.parser.parse(text)
        return self.transform.transform(tree)

    def write(self, data: dict) -> str:
        lines = []
        for key, val in data.items():
            if isinstance(val, list) and (len(val) == 0 or isinstance(val[0], list)):
                lines.append(f"\nbegin {key}")
                for row in val:
                    if not isinstance(row, list):
                        raise ValueError(f"Block '{key}' must be a list of lists; got row={row!r}")
                    lines.append(" ".join(_fmt_number(x) for x in row))
                lines.append("end")
            elif val is None:
                lines.append(f"{key}")
            elif isinstance(val, (int, float)):
                lines.append(f"{key} {_fmt_number(val)}")
            else:
                lines.append(f"{key} {val}")

        out = "\n".join(lines)
        if not out.endswith("\n"):
            out += "\n"
        return out
#endregion
