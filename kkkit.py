#
# from omegagene

import sys
import struct as st

def eliminate_comment(line, comment_char=[";","#"]):
    i = 0
    comment_index = -1
    for cc in comment_char:
        try:
            i = line.index(cc)
            if i >= 0 and (comment_index == -1 or i < comment_index):
                comment_index = i
        except ValueError:
            pass
    if comment_index >= 0:
        line = line[0:comment_index]
    return line

class FileIO(object):
    def __init__(self, fn):
        self.fn = fn
        self.f = None
    def open(self,mode):
        self.f = open(self.fn, mode)
        return self.f
    def close(self):
        return self.f.close()

class FileI(FileIO):
    def __init__(self, fn):
        super(FileI, self).__init__(fn)
    def open(self):
        return super(FileI,self).open("r")
    def readline_comment(self):
        try:
            line = self.f.readline()
            if line=="":  return None
        except: return None
        i = len(line)
        try:
            i = line.index(';')
        except ValueError:
            i = len(line) 
        line = line[0:i] + '\n'
        if len(line.strip()) == 0:
            return self.readline_comment()
        return line
    def read_line_com(self, marks=[";","#"]):
        line = self.f.readline();
        if not line: return None
        line = eliminate_comment(line)
        ##print "DD: " + line
        if len(line)==0: line = " "
        return line

class FileO(FileIO):
    def __init__(self, fn):
        super(FileO, self).__init__(fn)
    def open(self):
        return super(FileO,self).open("w")

class FileBI(FileIO):
    def __init__(self, fn):
        self.endian = ""
        super(FileBI, self).__init__(fn)
    def open(self):
        return super(FileBI,self).open("rb")
    def read_value(self, read_type, n_values=1):
        size = 4
        if read_type == "d": size=8
        read_st = self.endian + str(n_values) + read_type
        buf = self.f.read(size*n_values)
        val = st.unpack(read_st, buf)
        return val
    def read_int(self, n_values=1):
        return self.read_value("i", n_values)
    def read_float(self, n_values=1):
        return self.read_value("f", n_values)
    def read_double(self, n_values=1):
        return self.read_value("d",n_values)
        
class FileBO(FileIO):
    def __init__(self, fn):
        self.endian = ""
        super(FileBO, self).__init__(fn)
    def open(self):
        return super(FileBO,self).open("wb")
    
    
def err(msg):
    print (msg)
    return sys.stderr.write(msg)
