#region modules
from lark import Lark, Transformer
import logging 
#endregion

#region variables
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
#endregion

#region functions
#endregion

#region classes
class NamelistTransform(Transformer):
    def __init__(self):
        self.dict: dict = {}

    def file(self, args):
        return self.dict

    def namelist(self, args):
        logger.debug(f'namelist: {args}')
        name = args[1].value
        namelist_dict = {}
        for pair in args[2:-1]:
            key = pair[0].strip(' ')
            value = pair[1].strip(" '")

            namelist_dict.update({key: value})
        self.dict.update({name: namelist_dict})

    def pair(self, args):
        return (args[0].value, args[2].value)

class NamelistGrammar:
    grammar: str = r'''
!start: namelist*    -> file
!namelist: "&" NAME pair* "/" -> namelist
!pair: NAME "=" /\S+/ -> pair

NAME: /[a-zA-Z_0-9!]+/

%import common.WS
%ignore WS
'''
    transform: NamelistTransform = NamelistTransform()

    def __init__(self):
        self.parser = Lark(self.grammar, parser='lalr')

    def read(self, text: str) -> dict:
        tree = self.parser.parse(text)
        return self.transform.transform(tree)

    def write(self, data: dict):
        output = ''

        for namelist_name, namelist_dict in data.items():
            output += f'&{namelist_name}\n'
            for key, value in namelist_dict.items():
                output += f'{key}={value}\n'
            output += '/\n\n'

        return output
    
#endregion