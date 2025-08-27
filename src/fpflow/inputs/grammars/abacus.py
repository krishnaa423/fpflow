#region modules
from lark import Transformer, Lark
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class AbacusInputTransform(Transformer):
    pass

class AbacusStruTransform(Transformer):
    pass

class AbacusKptTransform(Transformer):
    pass

class AbacusInputGrammar:
    grammar = r'''
%import common.NEWLINE
%import common.WS_INLINE
%ignore WS_INLINE

NUMBER: /[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?/
NAME: /[A-Za-z_][A-Za-z0-9_\.]*/
VALUE: /.*/

start: "INPUT_PARAMETERS" NEWLINE pair+ NEWLINE*

pair: NAME VALUE NEWLINE
'''
    transform = AbacusInputTransform()

    def __init__(self):
        self.parser = Lark(self.grammar, parser='lalr')

    def read(self, text: str) -> dict:
        tree = self.parser.parse(text)
        return self.transform.transform(tree)

    def write(self, data: dict) -> str:
        pass

class AbacusStruGrammar:
    grammar = r'''
%import common.NEWLINE
%import common.WS_INLINE
%ignore WS_INLINE

NUMBER: /[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?/
SYMBOL: /[A-Z][a-z]?/
NAME: /[A-Za-z_][A-Za-z0-9_\.]*/

start: NEWLINE* items NEWLINE*

items: item+

item: atomic_species | numerical_orbital | lattice_constant | lattice_vectors | atomic_positions

atomic_species: "ATOMIC_SPECIES" NEWLINE species_line+ NEWLINE
species_line: SYMBOL NUMBER NAME NEWLINE

numerical_orbital: "NUMERICAL_ORBITAL" NEWLINE orbital_line+ NEWLINE
orbital_line: NAME NEWLINE

lattice_constant: "LATTICE_CONSTANT" NEWLINE NUMBER NEWLINE

lattice_vectors: "LATTICE_VECTORS" NEWLINE vectors_line+ NEWLINE
vectors_line: NUMBER NUMBER NUMBER NEWLINE

atomic_positions: "ATOMIC_POSITIONS" NEWLINE "Direct" NEWLINE position_entry+ NEWLINE
position_entry: SYMBOL NEWLINE NUMBER NEWLINE NUMBER NEWLINE position_vector+
position_vector: NUMBER NUMBER NUMBER NUMBER NUMBER NUMBER NEWLINE
'''
    transform = AbacusStruTransform()

    def __init__(self):
        self.parser = Lark(self.grammar, parser='lalr')

    def read(self, text: str) -> dict:
        tree = self.parser.parse(text)
        return self.transform.transform(tree)

    def write(self, data: dict) -> str:
        pass

class AbacusKptGrammar:
    grammar = r'''
%import common.NEWLINE
%import common.WS_INLINE
%ignore WS_INLINE

NUMBER: /[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?/
TYPE: "Line" | "Direct" | "Gamma"

start: "K_POINTS" NEWLINE NUMBER NEWLINE TYPE NEWLINE (gamma_line | kpt_line)+ NEWLINE*

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
        pass
#endregion