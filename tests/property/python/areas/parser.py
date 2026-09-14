"""Direct Sage parser/tokenizer comparisons used by polynomial string construction."""
import json
from sage.all import Integer
from sage.misc.parser import Parser, Tokenizer, LookupNameMaker, token_to_str

def comparison(run):
    def normalize(x):
        if isinstance(x,(list,tuple)): return [normalize(y) for y in x]
        if isinstance(x,bool) or x is None: return x
        return str(x)
    try: return json.dumps({'value':normalize(run())},separators=(',',':'),ensure_ascii=False)
    except Exception as e: return json.dumps({'error':type(e).__name__,'message':str(e)},separators=(',',':'),ensure_ascii=False)

def tokens(codes,mode,position):
    t=Tokenizer(''.join(chr(c) for c in codes))
    if mode==0: return [t.test(),t.test()]
    if mode==1:
        before=t.peek(); first=t.next(); last=t.last(); value=t.last_token_string(); after=t.peek()
        return [token_to_str(before),token_to_str(first),token_to_str(last),value,token_to_str(after),t.test()]
    if mode==2:
        t.test(); t.reset(position); return t.test()
    if mode==3:
        t.next(); t.backtrack(); t.backtrack(); return t.test()
    first=t.next(); result=t.backtrack(); again=t.next()
    return [token_to_str(first),result,token_to_str(again),t.last_token_string(),t.test()]

class Node:
    def __init__(self,s): self.s=s
    def __str__(self): return self.s
    def __neg__(self): return Node('(- '+str(self)+')')
    def binary(self,op,other): return Node('('+op+' '+str(self)+' '+str(other)+')')
for name,op in [('add','+'),('sub','-'),('mul','*'),('truediv','/'),('pow','^'),('eq','='),('ne','NOT_EQ'),('lt','<'),('le','LESS_EQ'),('gt','>'),('ge','GREATER_EQ')]:
    setattr(Node,'__'+name+'__',lambda self,other,op=op:self.binary(op,other))

def ast(codes,mode,implicit):
    source=''.join(chr(c) for c in codes)
    p=Parser(lambda s:Node('I'+str(Integer(s))),lambda s:Node('F'+s),lambda s:Node('V'+s),implicit_multiplication=bool(implicit))
    if mode==8: return str(p._variable_constructor()(source))
    if mode==9: return p._callable_constructor()(source)
    if mode==0: return str(p.parse(source))
    if mode==1: return str(p.parse_expression(source))
    t=Tokenizer(source)
    name=['p_eqn','p_expr','p_term','p_factor','p_power','p_atom'][mode-2]
    result=getattr(p,name)(t)
    return [str(result),t.test()]

def lookup(codes,fallback):
    name=''.join(chr(c) for c in codes)
    maker=LookupNameMaker({'a':'old','constructor':'owned'},(lambda n:'fallback:'+n) if fallback else None)
    first=maker(name)
    maker.set_names({'a':'new','__proto__':'ordinary'})
    return [first,maker(name)]

FUNCTIONS={
 'parser_tokens':lambda *args:comparison(lambda:tokens(*args)),
 'parser_ast':lambda *args:comparison(lambda:ast(*args)),
 'parser_lookup':lambda *args:comparison(lambda:lookup(*args)),
 'parser_token_name':lambda n:comparison(lambda:token_to_str(n)),
}

def constructed_names(mode):
    p=Parser(Node,Node,{'a':Node('A')},lambda name:lambda arg:Node(name+':'+str(arg)))
    return str(p.parse('a+a')) if mode==0 else str(p._callable_constructor()('foo')(Node('argument')))
FUNCTIONS['parser_constructed_names'] = lambda mode:comparison(lambda:constructed_names(mode))
