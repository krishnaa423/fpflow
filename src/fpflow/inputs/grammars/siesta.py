#region modules
from lark import Transformer, Lark
import re 
#endregion

#region variables
_TRUE = {"true",".true.","t","y","yes"}
_FALSE= {"false",".false.","f","n","no"}
#endregion

#region functions
def _coerce(tok: str):
    s = tok.strip()
    ls = s.lower()
    if ls in _TRUE:  return True
    if ls in _FALSE: return False
    # normalize Fortran D exponents
    s2 = re.sub(r'[dD]([+-]?\d+)$', r'e\1', s)
    if re.fullmatch(r'[+-]?\d+', s2): return int(s2)
    if re.fullmatch(r'[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:e[+-]?\d+)?', s2):
        try: return float(s2)
        except ValueError: pass
    return s
#endregion

#region classes
class SiestaTransform(Transformer):
    # unwrap the intermediary rule so start() gets (key, value) tuples
    def item(self, children):
        return children[0]

    # collect (k, v) pairs into a dict
    def start(self, items):
        return {k: v for (k, v) in items}

    # map: NAME VALUE
    def map(self, items):
        key = str(items[0])
        rhs = str(items[1])          # your VALUE token is the rest of the line
        toks = rhs.split()           # use shlex.split(rhs) if you need quoted strings
        vals = [_coerce(t) for t in toks]
        return (key, vals[0] if len(vals) == 1 else vals)

    # block: "%block" NAME NEWLINE block_line+ "%endblock" NAME NEWLINE
    def block(self, items):
        name = None
        rows = []
        for obj in items:
            # first NAME we see is the block name
            if name is None and getattr(obj, "type", None) == "NAME":
                name = str(obj)
            # each block_line already returns a list
            elif isinstance(obj, list):
                rows.append(obj)
        return (name, rows)

    # block_line: (NUMBER | SYMBOL)+ NEWLINE
    def block_line(self, items):
        # items are NUMBER/SYMBOL tokens (NEWLINE doesn’t arrive here)
        return [_coerce(str(t)) for t in items]

class SiestaGrammar:
    grammar: str = r"""
%import common.NEWLINE
%import common.WS_INLINE
%import common.ESCAPED_STRING -> STRING
%ignore WS_INLINE

COMMENT: /#[^\n]*/
%ignore COMMENT

NAME: /[A-Za-z][A-Za-z0-9_\.]*/
VALUE: /[^\n]+/                // rest of line, not including newline
SYMBOL: /[A-Z][a-z]?/
NUMBER: /[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?/
WORD: /[^\s#%][^\s#]*/         // generic token for units/labels (Ry, Ang, DZP, gamma, etc.)

// allow blank lines anywhere and EOF w/o final newline
_NL: (NEWLINE)+
start: _NL* (item (_NL+ item)*) _NL*

item: map | block

map: NAME VALUE                 // trailing newline handled by start separators

BLOCK_KW: /%[bB][lL][oO][cC][kK]/
ENDBLOCK_KW: /%[eE][nN][dD][bB][lL][oO][cC][kK]/

block: BLOCK_KW NAME _NL block_line+ ENDBLOCK_KW NAME?   // allow optional repeat name
block_line: (NUMBER | SYMBOL | WORD | STRING)+ _NL
"""

    transform: SiestaTransform = SiestaTransform()

    def __init__(self):
        self.parser = Lark(self.grammar, parser='lalr')

    def read(self, text: str) -> dict:
        tree = self.parser.parse(text)
        # return tree 
        return self.transform.transform(tree)
    
    def write(self, data: dict) -> str:
        def fmt(x):
            if isinstance(x, bool):
                return "T" if x else "F"
            if isinstance(x, int):
                return str(x)
            if isinstance(x, float):
                return f"{x:.10g}"
            return str(x)

        lines = []
        for key, val in data.items():
            # block: list of lists
            if isinstance(val, list) and (len(val) == 0 or isinstance(val[0], (list, tuple))):
                lines.append(f"\n%block {key}")
                for row in val:
                    if not isinstance(row, (list, tuple)):
                        raise ValueError(f"Block '{key}' must be a list of lists; got {row!r}")
                    lines.append(" ".join(fmt(x) for x in row))
                lines.append(f"%endblock {key}\n")
            # flag: None
            elif val is None:
                lines.append(key)
            # multi-value map line
            elif isinstance(val, (list, tuple)):
                lines.append(f"{key} " + " ".join(fmt(x) for x in val))
            # single-value map line
            else:
                lines.append(f"{key} {fmt(val)}")

        out = "\n".join(lines)
        if not out.endswith("\n"):
            out += "\n"
        return out

#endregion