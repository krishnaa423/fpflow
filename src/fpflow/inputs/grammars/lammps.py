#region modules
from lark import Lark, Transformer, Token
from fpflow.io.logging import get_logger
#endregion

#region variables
logger = get_logger()
#endregion

#region functions
#endregion

#region classes
class LammpsTransform(Transformer):
    def __init__(self):
        super().__init__(visit_tokens=True)
        self.dict: dict = {}

    def NAME(self, args):
        return args.value
    
    def VALUE(self, args):
        return args.value
    
    def pair(self, args):
        return args[:-1]

    def start(self, args):
        args = list(filter(lambda item: not isinstance(item, Token), args))

        for pair in args:
            key = pair[0]
            value = pair[1:]
            
            self.dict.update({key: value})

        return self.dict

class LammpsGrammar:
    grammar = r'''
%import common.NEWLINE
%import common.CNAME -> NAME
%import common.WS_INLINE
%ignore WS_INLINE

VALUE: /\S+/

start: NEWLINE* pair+

pair: NAME VALUE+ NEWLINE
'''
    transform: LammpsTransform = LammpsTransform()
    
    def __init__(self):
        self.parser: Lark = Lark(grammar=self.grammar, parser='lalr')

    def read(self, text: str) -> dict:
        tree = self.parser.parse(text)
        return self.transform.transform(tree)

    def write(self, data: dict) -> str:
        output = ''

        for key, value in data.items():
            output += f'{key} '
            if value is not None:
                for list_item in value:
                    output += f'{list_item} '

            output += f'\n'

        return output 

#endregion