"""Disable network connections in child calculations; installation is separate."""
import os,socket
if os.environ.get('V52_OFFLINE') == '1':
    def _blocked(*args,**kwargs):
        raise RuntimeError('Network connections are disabled for V52-R1 calculations')
    socket.create_connection=_blocked
    socket.socket.connect=_blocked
    socket.socket.connect_ex=_blocked
