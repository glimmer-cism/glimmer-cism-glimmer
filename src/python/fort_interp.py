import sys
#############################################
# having parsed the file with fort_parser.py,
# this code returns a structure describing the
# contents in terms of modules, programs, etc
#############################################

Disallowed='DisallowedFortTest'
FortSyntaxError='Fortran Syntax Error'
FortParseInternalError='Parser Internal Error'

printcoms=False

simple_types=('integer','real','complex','logical')

def print_comments(com):
    if printcoms: print 'Comment:',com

def dimsfind(l):
    out=[]
    for a in l:
        tmpdim=[]
        if a==',':
            out.append(tmpdim)
            tmpdim=[]
        else:
            tmpdim.append(a)
    out.append(tmpdim)
    return out

class FortTemplate:
    # Base class for others
    def interpmod(self,lines,next):
        if next>=len(lines):
            return next,False
        if lines[next].statements==[]:
            return next,False
        if lines[next].ismodule():
            out=Module()
            contains=False
            out.name=lines[next].statements[1]
            next=next+1
            while next<len(lines):
                if lines[next].isendmodule(out.name):
                    next=next+1
                    break
                if not contains:
                    next,found,decls=out.interpdecl(lines,next)
                    if found:
                        out.modvars=out.modvars+decls
                        continue
                    next,found=out.interptype(lines,next)
                    if found: continue
                    next,contains=out.interpcontains(lines,next)
                    if contains: continue
                else:
                    next,found=out.interpsub(lines,next)
                    if found: continue
                    next,found=out.interpfunc(lines,next)
                    if found: continue
                next,found=out.interpempty(lines,next)
                if found: continue
                next,found=out.interpunsupported(lines,next)
                if found: continue
            self.modules.append(out)
            return next,True
        return next,False

    ####################################
    
    def interpprog(self,lines,next):
        if next>=len(lines):
            return next,False
        if lines[next].statements==[]:
            return next,False
        if lines[next].isprogram():
            out=Program()
            contains=False
            out.name=lines[next].progname()
            next=next+1
            while next<len(lines):
                if lines[next].isendprogram(out.name):
                    next=next+1
                    break
                if not contains:
                    next,found,decls=out.interpdecl(lines,next)
                    if found:
                        out.decls=out.decls+decls
                        continue
                    next,found=out.interptype(lines,next)
                    if found: continue
                    next,contains=out.interpcontains(lines,next)
                    if contains: continue
                else:
                    next,found=out.interpsub(lines,next)
                    if found: continue
                    next,found=out.interpfunc(lines,next)
                    if found: continue
                next,found=out.interpempty(lines,next)
                if found: continue
                next,found=out.interpunsupported(lines,next)
                if found: continue
            self.progs.append(out)
            return next,True
        return next,False

    ####################################
    
    def interpsub(self,lines,next):
        if next>=len(lines):
            return next,False
        if lines[next].statements==[]:
            return next,False
        if lines[next].issub():
            out=Sub()
            contains=False
            out.name=lines[next].subname()
            argnames=lines[next].subargnames()
            next=next+1
            while next<len(lines):
                if lines[next].isendsub(out.name):
                    next=next+1
                    break
                if not contains:
                    next,found,decls=out.interpdecl(lines,next)
                    if found:
                        for a in decls:
                            if a.name in argnames:
                                out.args.append(a)
                            else:
                                out.intvars.append(a)
                        continue
                    next,found=out.interptype(lines,next)
                    if found: continue
                    next,contains=out.interpcontains(lines,next)
                    if contains: continue
                else:
                    next,found=out.interpsub(lines,next)
                    if found: continue
                    next,found=out.interpfunc(lines,next)
                    if found: continue
                next,found=out.interpempty(lines,next)
                if found: continue
                next,found=out.interpunsupported(lines,next)
                if found: continue
            self.subs.append(out)
            return next,True
        return next,False

    ####################################
    
    def interpfunc(self,lines,next):
        if next>=len(lines):
            return next,False
        if lines[next].statements==[]:
            return next,False
        return next,False

    ####################################
    
    def interpdecl(self,lines,next):
        out=[]
        if next>=len(lines):
            return next,False,out
        if lines[next].statements==[]:
            return next,False,out
        if lines[next].isdecl():
            tmpdecl=lines[next].getdecl()
            for a in tmpdecl['varnames']:
                out.append(Decl())
                out[-1].name=a
                out[-1].type=tmpdecl['type']
                if tmpdecl.has_key('kind'):
                    out[-1].kind=tmpdecl['kind']
                if tmpdecl.has_key('len') and out[-1].type=='character':
                    out[-1].len=tmpdecl['len']
                else:
                    out[-1].len='1'
                if tmpdecl.has_key('attributes'):
                    out[-1].allocatable=tmpdecl['attributes'].get('allocatable',False)
                    out[-1].intent=tmpdecl['attributes'].get('intent',[''])[0]
                    out[-1].parameter=tmpdecl['attributes'].get('parameter',False)
                    out[-1].pointer=tmpdecl['attributes'].get('pointer',False)
                    out[-1].public=tmpdecl['attributes'].get('public',True)
                    out[-1].public=(not tmpdecl['attributes'].get('private',False))
                    if 'dimension' in tmpdecl['attributes']:
                        out[-1].array=True
                        out[-1].size=dimsfind(tmpdecl['attributes']['dimension'])
                out[-1].comments=lines[next].comment
            next=next+1
            return next,True,out
        return next,False,out

    ####################################

    def interpcontains(self,lines,next):
        if next>=len(lines):
            return next,False
        if lines[next].statements==[]:
            return next,False
        if lines[next].iscontains():
            return next+1,True
        return next,False

    ####################################
     
    def interptype(self,lines,next):
        if next>=len(lines):
            return next,False
        if lines[next].statements==[]:
            return next,False
        if lines[next].istypedef():
            out=Type()
            out.name=lines[next].typename()
            while next<len(lines):
                if lines[next].isendtypedef(out.name):
                    next=next+1
                    break
                next,found,decls=out.interpdecl(lines,next)
                if found:
                    out.components=out.components+decls
                    continue
                next,found=out.interpempty(lines,next)
                if found: continue
                next,found=out.interpunsupported(lines,next)
                if found: continue
            self.typedefs.append(out)
            return next,True
        return next,False

    ####################################
    
    def interpempty(self,lines,next):
        if next>=len(lines):
            return next,False
        if lines[next].hascomment():
            self.comments.append(lines[next].comment)
            return next+1,True
        else:
            return next,False

    ####################################

    def interpunsupported(self,lines,next):
        if next>=len(lines):
            return next,False
        next=next+1
        return next,True

class FortFile(FortTemplate):
    def __init__(self):
        self.modules=[]
        self.progs=[]
        self.subs=[]
        self.funcs=[]
        self.comments=[]
    def interpdecl(self,lines,next):
        raise Disallowed
    def interptype(self,lines,next):
        raise Disallowed
    def interpcontains(self,lines,next):
        raise Disallowed
    def Print(self):
        for a in self.comments:
            print_comments(a)
        for a in self.modules:
            a.Print()
        for a in self.progs:
            a.Print()
        for a in self.subs:
            a.Print()

class Program(FortTemplate):
    def __init__(self):
        self.name=''
        self.subs=[]
        self.funcs=[]
        self.decls=[]
        self.types=[]
        self.comments=[]
    def interpmod(self,lines,next):
        raise Disallowed
    def interpprog(self,lines,next):
        raise Disallowed
    def Print(self):
        print 'Program',self.name
        for a in self.comments:
            print_comments(a)
        for a in self.subs:
            a.Print()

class Module(FortTemplate):
    def __init__(self):
        self.name=''
        self.modvars = []
        self.typedefs = []
        self.subs = []
        self.funcs = []
        self.comments=[]
    def interpmod(self,lines,next):
        raise Disallowed
    def interpprog(self,lines,next):
        raise Disallowed
    def Print(self):
        print 'Module',self.name
        for a in self.comments:
            print_comments(a)
        for a in self.modvars:
            a.Print()
        for a in self.typedefs:
            a.Print()
        for a in self.subs:
            a.Print()

class Type(FortTemplate):
    def __init__(self):
        self.name=''
        self.components = []
        self.comments=[]
    def interpmod(self,lines,next):
        raise Disallowed
    def interpprog(self,lines,next):
        raise Disallowed
    def interpsub(self,lines,next):
        raise Disallowed
    def interpfunc(self,lines,next):
        raise Disallowed
    def interptype(self,lines,next):
        raise Disallowed
    def interpcontains(self,lines,next):
        raise Disallowed
    def Print(self):
        print 'Type definition',self.name
        for a in self.comments:
            print_comments(a)
        for a in self.components:
            a.Print()

class Sub(FortTemplate):
    def __init__(self):
        self.name=''
        self.args=[]
        self.intvars=[]
        self.typedefs = []
        self.subs=[]
        self.funcs=[]
        self.comments=[]
    def interpmod(self,lines,next):
        raise Disallowed
    def interpprog(self,lines,next):
        raise Disallowed
    def Print(self):
        print 'Subroutine',self.name,'('+','.join(self.argnames())+')'
        for a in self.comments:
            print_comments(a)
        for a in self.subs:
            a.Print()
    def argnames(self):
        out=[]
        for a in self.args:
            out.append(a.name)
        return out

##############################

class Decl:
    def __init__(self):
        self.name=''
        self.type=''
        self.pointer=False
        self.array=False
        self.optional=False
        self.parameter=False
        self.allocatable=False
        self.public=True
        self.save=True
        self.size=[]
        self.kind=''
        self.len=''
        self.intent=''
        self.comments=[]
    def Print(self):
        out=''
        if self.type in simple_types:
            out=out+self.type+{'':''}.get(self.kind,'('+self.kind+')')
        elif self.type=='character':
            out=self.type+'('+self.len+{'':''}.get(self.kind,',kind='+self.kind)+')'
        else:
            out='type('+self.type+')'
        if self.pointer:   out=out+',pointer'
        if self.parameter: out=out+',parameter'
        if self.optional:  out=out+',optional'
        if self.allocatable:  out=out+',allocatable'
        if (not self.public): out=out+',private'
        if self.array:     out=out+',array(rank='+str(len(self.size))+')'
        if self.intent!='': out=out+',intent('+self.intent+')'
        out=out+' :: '+self.name
        print out

##############################
# Main call
##############################

def InterpFortFile(lines):

    out=FortFile()
    next=0   # Next line to be interpreted
    while next<len(lines):
        next,found=out.interpmod(lines,next)
        if found: continue
        next,found=out.interpprog(lines,next)
        if found: continue
        next,found=out.interpsub(lines,next)
        if found: continue
        next,found=out.interpfunc(lines,next)
        if found: continue
        next,found=out.interpempty(lines,next)
        if found: continue
        next,found=out.interpunsupported(lines,next)
    return out
