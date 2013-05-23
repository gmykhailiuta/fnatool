#!/usr/bin/python
#from process import process_signal
from gi.repository import Gtk

class MyWindow(Gtk.Window):
    def __init__(self):
        Gtk.Window.__init__(self, title="Hello World")

        self.box = Gtk.Box(spacing=6)
        self.add(self.box)
        _pack = self.box.pack_start

        _pack(Gtk.Label(label='Record'), True, True, 0)
        self.ent_record = Gtk.Entry(text="")
        #self.ent_record.set_text("Hello World")
        _pack(self.ent_record, True, True, 0)

        #self.lb_window = Gtk.Label(label='Window')
        _pack(Gtk.Label(label='Window'), True, True, 0)
        self.sbt_window = Gtk.SpinButton()
        _pack(self.sbt_window, True, True, 0)
        
        self.btn_process = Gtk.Button(label="Process")
        self.btn_process.connect("clicked", self.on_btn_process_clicked)
        _pack(self.btn_process, True, True, 0)


        self.btn_close = Gtk.Button(label="Close")
        self.btn_close.connect("clicked", self.on_btn_close_clicked)
        self.box.pack_start(self.btn_close, True, True, 0)

    def on_btn_process_clicked(self, widget):
        print "Process"
        #process_signal(self.ent_record.get_text())

    def on_btn_close_clicked(self, widget):
        Gtk.main_quit()

win = MyWindow()
win.connect("delete-event", Gtk.main_quit)
win.show_all()
Gtk.main()
