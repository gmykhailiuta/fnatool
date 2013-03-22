#!/usr/bin/env python
import os

fname = "beta.htm"
fh = open(fname, "w")
fh.write("<html><head></head><body>")
fh.write("<table>")
for filename in sorted(os.listdir('.')):
    path = os.path.splitext(filename)
    if path[1] == '.hea':
        fh.write("<tr><td>")
        for c in range(2):
            fh.write('<img width="400" height="400" src=%s_c%s_std.png /><img width="400" height="400" src=%s_c%s_sdn.png /><br/>' % (path[0],c,path[0],c))
        fh.write("</td><td valign=\"top\">")
        hh = open(filename, "r")
        header = hh.read()
        header = header.replace("\n", "<br/>")
        header = header.replace(" ", "&nbsp;")
        header = header.replace("#", "")
        hh.close()
        fh.write(header)
        fh.write("</td></tr><tr height=2 bgcolor=black><td></td><td></td></tr>")
fh.write("</table></body>")

                
