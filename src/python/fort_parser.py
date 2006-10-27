#!/usr/bin/python
# An attempt to create a general fortran parser written in Python
# by Ian Rutt

from fort_interp import *

name_chars='abcdefghijklmnopqrstuvwxyz_0123456789'
number_chars='0123456789'
exponents='eEdD'
seven_char_tokens=['.false.']
six_char_tokens=['.true.']
five_char_tokens=['.not.','.and.']
four_char_tokens=['.le.','.gt.','.lt.','.ge.','.eq.','.ne.','.or.']
three_char_tokens=[]
two_char_tokens=['::','<=','>=','==','=>','/)','(/','//','**','/=','""',"''"]
one_char_tokens=['-','_','+','=','/','*',' ','>','<','(',')',',',':','%']
types=['integer','real','complex','character','logical']
function_words=['function','recursive','pure','elemental']
attribute_names=['dimension','optional','pointer','allocatable','intent','parameter','save']
allowed_intents=['in','out','inout']
suspect_number_endings=['.e','.d']

##### Syntax checking

def valid_name(name):
    return (name[0].isalpha() and name.rstrip(name_chars)=='')

def names_only(l):
    out=[]
    for a in l:
        if valid_name(str(a)):
            out.append(a)
    return out

def handle_var_params(l):
    if len(l)==1:
        return l[0]
    elif len(l)==3 and l[0:3]==['kind','=']:
        return l[2]
    else:
        raise FortSyntaxError

def get_intent_contents(l):
    if l[0]!='(' or l[2]!=')':
        raise FortSyntaxError
    elif l[1] not in allowed_intents:
        raise FortSyntaxError
    else:
        return l[1],3

def get_bracket_contents(l):
    if l[0]=='(':
        return l[1:l.index(')')],l.index(')')+1
    else:
        raise FortSyntaxError

#############################################################
# Processing declarations - functions
#############################################################

def typereal(l):
    decldata={'()':genkind,
              ',':attribs,
              '::':varnames}.get(str(l[0]),varnames)(l)
    decldata.update({'type':'real'})
    return decldata

def typeint(l):
    decldata={'()':genkind,
              ',':attribs,
              '::':varnames}.get(str(l[0]),varnames)(l)
    decldata.update({'type':'integer'})
    return decldata

def typechar(l):
    decldata={'()':charkind,
              ',':attribs,
              '::':varnames}.get(str(l[0]),varnames)(l)
    decldata.update({'type':'character'})
    return decldata

def typecomplex(l):
    decldata={'()':genkind,
              ',':attribs,
              '::':varnames}.get(str(l[0]),varnames)(l)
    decldata.update({'type':'complex'})
    return decldata

def typelogical(l):
    decldata={'()':genkind,
              ',':attribs,
              '::':varnames}.get(str(l[0]),varnames)(l)
    decldata.update({'type':'logical'})
    return decldata


def typetype(l):
    typename=l[0].contents[0]
    decldata={',':attribs,
              '::':varnames}.get(str(l[1]),varnames)(l[1:])
    decldata.update({'type':typename})
    return decldata

def genkind(l):
    if len(l[0].contents)==3:
        kind=l[0].contents[2]
    else:
        kind=l[0].contents[0]
    l=l[1:]
    decldata={',':attribs,
              '::':varnames}.get(str(l[0]),varnames)(l)
    decldata.update({'kind':kind})
    return decldata

def charkind(l):
    kindout={}
    c=l[0].contents
    try:
        c=[c[:c.index(',')],c[c.index(',')+1:]]
    except ValueError:
        c=[c]
    except:
        raise
    if len(c)>2:
        raise FortSyntaxError,c
    for a in c:
        if len(a)==1:
            kindout.update({'len':a[0]})
        elif len(a)==3:
            kindout.update({a[0]:a[2]})
    l=l[1:]
    decldata={',':attribs,
              '::':varnames}.get(str(l[0]),varnames)(l)
    decldata.update(kindout)
    return decldata

def attribs(l):
    att={}
    while 1:
        if l[0]==',':
            l=l[1:]
            continue
        elif l[0]=='optional':
            att.update({'optional':True})
            l=l[1:]
            continue
        elif l[0]=='pointer':
            att.update({'pointer':True})
            l=l[1:]
            continue
        elif l[0]=='intent':
            att.update({'intent':l[1].contents})
            l=l[2:]
            continue
        elif l[0]=='allocatable':
            att.update({'allocatable':True})
            l=l[1:]
            continue
        elif l[0]=='dimension':
            att.update({'dimension':l[1].contents})
            l=l[2:]
            continue
        elif l[0]=='parameter':
            att.update({'parameter':True})
            l=l[1:]
            continue
        elif l[0]=='private':
            att.update({'private':True})
            l=l[1:]
            continue
        elif l[0]=='public':
            att.update({'public':True})
            l=l[1:]
            continue
        elif l[0]=='save':
            att.update({'save':True})
            l=l[1:]
            continue
        elif l[0]=='':
            l=l[1:]
            continue
        elif l[0]=='::':
            break
        elif valid_name(str(l[0])):
            break
        else:
            SynErr(l)
        
    vn=varnames(l)
    vn.update({'attributes':att})
    return vn

def SynErr(line):
    raise FortSyntaxError,line

def varnames(l):
    if l[0]=='::': l=l[1:]
    vn=[]
    i=0
    assign=False
    while i<len(l):
        if (not assign) and valid_name(str(l[i])):
            vn.append(l[i])
            i=i+1
        elif l[i]==',':
            assign=False
            i=i+1
        else:
            assign=True
            i=i+1
    
    return {'varnames':vn}

##### For removing delimited sections

def remove_delimited(line,d1,d2,DoubleCode=False):

    # Set DoubleCode to allow double delimiters to be retained

    bk=0
    temp_str=''
    undel_str=''
    delimited=[]

    i=0
    while i<len(line):
        if bk==0:
            # If we're in an undelimited section
            if line[i]==d1:
                bk=1
                test=1
            undel_str=undel_str+line[i]
            i=i+1
            continue
        elif bk==1:
            # If we're in a delimited section
            if line[i]!=d1 and line[i]!=d2:
                # Still in the section
                temp_str=temp_str+line[i]
                i=i+1
                continue
            elif DoubleCode:
                # If we're testing for doublecount
                try:
                    tmp=line[i+1]
                except IndexError:
                    tmp=None
                except:
                    raise
                if line[i]==tmp and (line[i]==d1 or line[i]==d2):
                    # Found double delimiter
                    temp_str=temp_str+line[i]+line[i+1]
                    i=i+2
                    continue
                else:
                    # Found single delimiter
                    if line[i]==d2:
                        delimited.append(temp_str[:])
                        bk=0
                        temp_str=''
                        undel_str=undel_str+line[i]
                        i=i+1
                        continue
                    else:
                        raise FortSyntaxError
            else:
                # Normal circumstances
                if line[i]==d2:
                    delimited.append(temp_str[:])
                    bk=0
                    temp_str=''
                    undel_str=undel_str+line[i]
                    i=i+1
                    continue
                else:
                    temp_str=temp_str+line[i]
                    i=i+1
                    continue

    return delimited,undel_str

#### For putting them back

def recover_delimited(line,d1,d2,delimited):

    if delimited==[]:
        return line,[]
    
    i=0
    while i<len(line):
        if line[i]==d1:
            line=line[0:i+1]+delimited[0]+line[i+1:]
            i=i+len(delimited[0])+1
            delimited=delimited[1:]
        if delimited==[]:
            break
        i=i+1

    return line,delimited

############################################################

class Brackets:

    def __init__(self,contents):
        self.isbrackets=True
        self.isarrayconstruct=False
        self.contents=contents

    def __str__(self):
        return '()'

class ArrayConstruct:

    def __init__(self,contents):
        self.isbrackets=False
        self.isarrayconstruct=True
        self.contents=contents
        
    def __str__(self):
        return '(//)'

def newbracket(line,i):

    if line[i]=='(' or line[i]=='(/':
        if line[i]=='(':  out=Brackets([])
        if line[i]=='(/': out=ArrayConstruct([])
        i=i+1
        while i<len(line):
            if line[i]=='(' or line[i]=='(/':
                nb,i=newbracket(line,i)
                out.contents.append(nb)
            elif line[i]==')' and out.isbrackets:
                i=i+1
                return out,i
            elif line[i]=='/)' and out.isarrayconstruct:
                i=i+1
                return out,i
            else:
                out.contents.append(line[i])
                i=i+1
    else:
        i=i+1

def Bracketize(line):
    result=[]
    i=0
    while i<len(line):
        if line[i]=='(' or line[i]=='(/':
            nb,i=newbracket(line,i)
            result.append(nb)
        else:
            result.append(line[i])
            i=i+1
    return result

############################################################

class FortLine:

    def __init__(self):
        self.comment=''
        self.statements=''
        self.label=''
        self.continues=False
        self.continuing=False
        self.line=''

    def Print(self):
        print self.label,self.statements,
        if self.comment.strip()!='':
            print '!',self.comment
        else:
            print
        
    def __add__(self,other):
        ret=FortLine()
        ret.comment=self.comment+' '+other.comment
        ret.statements=self.statements+other.statements
        ret.label=self.label
        ret.continues=other.continues
        ret.continuing=self.continuing
        ret.line=self.line+other.line
        return ret

    def split(self):
        outlist=[]
        temp_state=self.statements
        dquote,temp_state=remove_delimited(temp_state,'"','"',True)
        squote,temp_state=remove_delimited(temp_state,"'","'",True)
        list_state=temp_state.split(';')
        for i in range(len(list_state)):
                list_state[i],dquote=recover_delimited(list_state[i],'"','"',dquote)
                list_state[i],squote=recover_delimited(list_state[i],"'","'",squote)
                out_temp=FortLine()
                out_temp.comment=self.comment
                out_temp.statements=list_state[i].strip()
                if i==0:
                    out_temp.label=self.label
                    out_temp.continuing=self.continuing
                    out_temp.line=self.line
                if i==len(list_state):
                    out_temp.continues=self.contines
                outlist.append(out_temp)

        return outlist        
    
    ### General stuff

    def labelcheck(self):
        temp_state=self.statements
        for i in range(len(temp_state)):
            if  not temp_state[0].isdigit(): break
            self.label=self.label+temp_state[0]
            temp_state=temp_state[1:].strip()
        self.statements=temp_state

    def empty(self):
        return (self.statements==[] and self.comment=='')

    def hascomment(self):
        return (self.comment.strip()!='')

    ### Modules

    def ismodule(self):
        return (self.statements[0]=='module')

    def isendmodule(self,name):
        if self.statements[0:2]==['end','module']:
            if self.statements[2]==name:
                return True
            else:
                raise FortSyntaxError

    ### Programs

    def isprogram(self):
        return (self.statements[0]=='program')

    def isendprogram(self,name):
        if self.statements[0:2]==['end','program']:
            if self.statements[2]==name:
                return True
            else:
                raise FortSyntaxError

    def progname(self):
        if self.isprogram():
            if len(self.statements)==2:
                return self.statements[1]
            else:
                return ''
        else:
            raise FortParseInternalError

    ### Subroutines

    def issub(self):
        return (self.statements[0]=='subroutine')

    def isendsub(self,name):
        if self.statements[0:2]==['end','subroutine']:
            if self.statements[2]==name:
                return True
            else:
                raise FortSyntaxError
    def subname(self):
        if self.issub():
            if len(self.statements)>=2 and valid_name(self.statements[1]):
                return self.statements[1]
            else:
                raise FortSyntaxError
        else:
            raise FortParseInternalError

    def subargnames(self):
        if self.issub():
            if len(self.statements)==2:
                return []
            elif str(self.statements[2])=='()':
                try:
                    return names_only(self.statements[2].contents)
                except ValueError:
                    raise FortSyntaxError
                except:
                    raise
            else:
                raise FortSyntaxError

    ### contains

    def iscontains(self):
        return (self.statements[0]=='contains')

    ### Declarations

    def isdecl(self):
        if self.statements[0] in types:
            if self.statements[1] in function_words:
                return False
            else:
                return True
        elif self.statements[0]=='type' and str(self.statements[1])=='()':
            if self.statements[2] in function_words:
                return False
            else:
                return True
        else:
            return False

    def getdecl(self):
        line=self.statements
        # Process each element in turn
        decldata={'character':typechar,
         'real':typereal,
         'integer':typeint,
         'complex':typecomplex,
         'logical':typelogical,
         'type':typetype}[line[0]](line[1:])
        return decldata

    ### Derived type definitions

    def istypedef(self):
        return (self.statements[0]=='type' and valid_name(str(self.statements[1])))

    def isendtypedef(self,name):
        if self.statements[0:2]==['end','type']:
            if self.statements[2]==name:
                return True
            else:
                raise FortSyntaxError
            
    def typename(self):
        if self.istypedef():
            if len(self.statements)>=2 and valid_name(str(self.statements[1])):
                return self.statements[1]
            else:
                raise FortSyntaxError
        else:
            raise FortParseInternalError

    
#############################################################

def Process(raw_lines,FixedFormat):

    # Empty list for processed lines
    processed_lines=[]

    # First pass
    for a in raw_lines:
        processed_lines.append(ProcessLine(a,FixedFormat))

    # Second pass
    processed_lines=ProcessLinesPart2(processed_lines)

    # Split into tokens
    for a in processed_lines:
        a.statements=Tokenize(a.statements)

    # Remove blank lines (i.e. no comments or statements)
    pl=[]
    for i in range(len(processed_lines)):
        if not processed_lines[i].empty():
            pl.append(processed_lines[i])

    # Process brackets
    for a in pl:
        a.statements=Bracketize(a.statements)
   
    return pl

#############################################################

def ProcessLine(raw,FixedFormat):

    out_line=FortLine()
    out_line.line=raw
    dquote=[]
    squote=[]
    raw=raw.rstrip()
    if raw=='':
        return out_line

    # Ignore preprocessor lines for now
    if raw[0]=='#':
        return out_line

    if(FixedFormat):
        for i in range(len(raw)):
            if i==0 and (raw[i]=='c' or raw[i]=='!'):
                out_line.comment=raw[1:].strip()
                return out_line
            if i>0 and i<5 and raw[i].isdigit():
                out_line.label+=raw[i]
            if i==5 and raw[i]!=' ':
                out_line.continuing=True
            if i==6:
                raw_state=raw[6:].strip()
    else:
        raw_state=raw

    # Remove quoted text so we can check for comments
    dquote,raw_state=remove_delimited(raw_state,'"','"',True)
    squote,raw_state=remove_delimited(raw_state,"'","'",True)

    if raw_state.find('!')!=-1:
        out_line.comment=raw_state.split('!')[-1]
        raw_state=raw_state.split('!')[0]

    # Put it back again
    raw_state,dquote=recover_delimited(raw_state,'"','"',dquote)
    raw_state,squote=recover_delimited(raw_state,"'","'",squote)
    out_line.comment,dquote=recover_delimited(out_line.comment,'"','"',dquote)
    out_line.comment,squote=recover_delimited(out_line.comment,"'","'",squote)

    # Check to see if this line continues
    if len(raw_state.strip())>0 and raw_state.strip()[-1]=='&':
        raw_state=raw_state.strip()[:-1].strip()
        out_line.continues=True

    out_line.statements=raw_state.strip()

    return out_line

#############################################################

def ProcessLinesPart2(linelist):

    # Joins continuation lines and their comments together

    i=0
    ii=len(linelist)
    outlist=[]

    while i<ii:
        thisline=linelist[i]
        i=i+1
        while i<ii and (thisline.continues or linelist[i].continuing):
            thisline=thisline+linelist[i]
            i=i+1
        outlist.append(thisline)

    outlist2=[]

    # Split up statements on the same line

    for a in outlist:
        outlist2=outlist2+a.split()

    # Find any remaining labels

    for a in outlist2:
        a=a.labelcheck()
        
    return outlist2

#############################################################

def checkident(state,list):
    s=state
    name=''
    if len(s)==0:
        return state,list
    if not s[0].isalpha():
        return state,list
    while len(s)!=0:
        if name_chars.find(s[0])!=-1:
            name=name+s[0]
            s=s[1:]
        else: break

    list.append(name)
    return s,list

#############################################################

def checknumber(state,list):
    s=state
    number=''
    exponent=False
    decimal=False
    if len(s)==0:
        return state,list
    if not (s[0].isdigit() or (s[0]=='.' and s[1].isdigit())):
        return state,list
    
    # Here we check for a number
    while len(s)!=0:
        if s[0].isdigit():
            number=number+s[0]
            s=s[1:]
        elif not decimal and s[0]=='.':
            number=number+s[0]
            s=s[1:]
            decimal=True
        elif exponents.find(s[0])!=-1:
            number=number+s[0]
            s=s[1:]
            exponent=True
            break
        else: break
    if exponent:
        if '+-'.find(s[0])!=-1:
            number=number+s[0]
            s=s[1:]
        while len(s)!=0:
            if s[0].isdigit():
                number=number+s[0]
                s=s[1:]
            else:
                break

    if number[-1]=='.' and len(s)!=0:
        if s[0].isalpha():
            s=number[-1]+s
            number=number[:-1]
    elif number[-2:] in suspect_number_endings:
        s=number[-2:]+s
        number=number[:-2]

    list.append(number)
    return s,list

#############################################################

def checktoken(state,list):
    s=state
    if len(s)==0:
        return state,list

    # Seven chars
    for a in seven_char_tokens:
        if s[0:7]==a:
            list.append(a)
            s=s[7:]
            return s,list
        
    # Six chars
    for a in six_char_tokens:
        if s[0:6]==a:
            list.append(a)
            s=s[6:]
            return s,list
        
    # Five chars
    for a in five_char_tokens:
        if s[0:5]==a:
            list.append(a)
            s=s[5:]
            return s,list
        
    # Four chars
    for a in four_char_tokens:
        if s[0:4]==a:
            list.append(a)
            s=s[4:]
            return s,list
        
    # Three chars
    for a in three_char_tokens:
        if s[0:3]==a:
            list.append(a)
            s=s[3:]
            return s,list
        
    # Two chars
    for a in two_char_tokens:
        if s[0:2]==a:
            if not s[0:3]=='(/=':
                # This is a bit of a fudge
                list.append(a)
                s=s[2:]
                return s,list
        
    # One char
    for a in one_char_tokens:
        if s[0:1]==a:
            s=s[1:]
            if a!=' ':
                list.append(a)
            return s,list

    return s,list

#############################################################

def Tokenize(state):
    # Splits a statement into identifiers and tokens

    s=state
    list=[]
    dquote,s=remove_delimited(s,'"','"',True)
    squote,s=remove_delimited(s,"'","'",True)

    while len(s)!=0:
        s,list=checkident(s,list)
        s,list=checknumber(s,list)
        s,list=checktoken(s,list)

    for i in range(len(list)):
        list[i],dquote=recover_delimited(list[i],'"','"',dquote)
        list[i],squote=recover_delimited(list[i],"'","'",squote)

    return list

#############################################################
# Top-level function call                                   #
#############################################################

def ParseFortran(filename):

    # Determine formatting
    extension=filename.split('.')[-1]
    try:
        FixedFormat={'f':True,'F':True,'f90':False,'F90':False}[extension]
    except:
        print 'Unrecognised filename extension'
        return 255

    # Open for input and read into list
    infile=open(filename,'r')
    raw_file=infile.readlines()
    infile.close()

    # Convert to lowercase
    # Replace tabs with spaces
    for i in range(len(raw_file)):
        raw_file[i]=raw_file[i].lower()
        raw_file[i]=raw_file[i].replace('\t','      ')

    # Process the lines in turn to sort out formatting
    processed_file=Process(raw_file,FixedFormat)
 
    return processed_file
