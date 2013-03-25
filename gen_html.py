#!/usr/bin/env python
import os

fname = "beta.htm"
fh = open(fname, "w")
fh.write("<html><head></head><body>")
fh.write("<table>")
for filename in sorted(os.listdir('.')):
    path = os.path.splitext(filename)
    if path[1] == '.hea':
        fh.write("<h2>%s</h2><table><tr><td>" % (path[0],))
        for c in range(2):
            fh.write("<h4>Channel %s</h4>" % (c+1,))
            fh.write('<table><tr><td><img width="400" height="400" src=%s_c%s_std.png /></td><td><img width="400" height="400" src=%s_c%s_sdn.png /></td></tr></table></div></p>' % (path[0],c,path[0],c))
        fh.write("</td><td valign=\"top\">")
        hh = open(filename, "r")
        header = hh.read()
        header = header.replace("\n", "<br/>")
        header = header.replace(" ", "&nbsp;")
        header = header.replace("#", "")
        hh.close()
        fh.write(header)
        fh.write("</td></tr></table><hr height=2/>")
fh.write("</body></html>")

                
