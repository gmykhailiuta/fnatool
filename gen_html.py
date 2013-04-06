#!/usr/bin/env python
import os

fname = "beta.htm"
fh = open(fname, "w")
fh.write("<html><head></head><body>")
for filename in sorted(os.listdir('.')):
	path = os.path.splitext(filename)
	if path[1] == '.hea':
		fh.write("<p style=\"page-break-before: always\"><h2>%s</h2><table><tr><td valign=\"top\">" % (path[0],))
		fh.write('<table><tr><td><img width="400" height="400" src=%s_std.png /></td></tr><tr><td><img width="400" height="400" src=%s_sdn.png /></td></tr></table>' % (path[0],path[0]))
		fh.write("</td><td valign=\"top\">")
		hh = open(filename, "r")
		header = hh.read()
		header = header.replace("\n", "<br/>")
		header = header.replace(" ", "&nbsp;")
		header = header.replace("#", "")
		hh.close()
		fh.write(header)
		fh.write("</td></tr></table><hr height=2 /></p>")
fh.write("</body></html>")

                
