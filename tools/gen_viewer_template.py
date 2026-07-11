#!/usr/bin/env python3
"""Regenerate viewer_template.h from tools/clannview.template.html.
Run from the repo root: python3 tools/gen_viewer_template.py"""
import os
here=os.path.dirname(os.path.abspath(__file__)); root=os.path.dirname(here)
html=open(os.path.join(here,'clannview.template.html')).read()
b='/*CLANN_DATA_BEGIN*/'; e='/*CLANN_DATA_END*/'
i=html.index(b)+len(b); j=html.index(e)
head=html[:i]; tail=html[j:]
def cstr(s):
    out=[]
    for line in s.splitlines(keepends=True):
        esc=line.replace('\\','\\\\').replace('"','\\"').replace('\n','\\n').replace('\t','\\t').replace('\r','')
        out.append('"'+esc+'"')
    return '\n'.join(out) if out else '""'
with open(os.path.join(root,'viewer_template.h'),'w') as f:
    f.write("/* Auto-generated from tools/clannview.template.html by tools/gen_viewer_template.py -- do not edit by hand. */\n")
    f.write("#ifndef VIEWER_TEMPLATE_H\n#define VIEWER_TEMPLATE_H\n\n")
    f.write("static const char *VIEWER_HTML_HEAD =\n"+cstr(head)+";\n\n")
    f.write("static const char *VIEWER_HTML_TAIL =\n"+cstr(tail)+";\n\n#endif\n")
print("regenerated viewer_template.h")
