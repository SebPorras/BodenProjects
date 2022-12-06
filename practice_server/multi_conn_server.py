###############################################################################
# Date: 2/12
# Author: Sebastian Porras   
# Aims: Create a simple multi-connection practice 
# server and client to practice socket programming.    
# Based off tutorial https://realpython.com/python-sockets/#reference
###############################################################################

import sys 
import socket 
import selectors 
import types 

#object that can handle multiple I/O completions 
#for more than one socket. 
sel = selectors.DefaultSelector()

#access the cli for host and port info 
host, port = (sys.argv[1], int(sys.argv[2]))

#create socket and configure with IPv4 and TCP protocol 
lsock = socket.socket((socket.AF_INET, socket.SOCK_STREAM))
lsock.bind((host, port))
lsock.listen()
print(f"Listening on {(host, port)}")
lsock.setblocking(False)

#register the socket to be monitored 
sel.register(lsock, selectors.EVENT_READ, data=None)


####Accept wrapper: Used when client socket has already ben accepted####

def accept_wrapper(sock):
    #conn (socket object) should be ready to read 
    conn, addr = sock.accept()
    print(f"Accepted connection from {addr}")
    conn.setblocking(False)
    data = types.SimpleNamespace(addr=addr, inb=b"", outb=b"")




#### EVENT LOOP ####

try:
    while True: 
        #blocks until sockets are ready for I/O
        events = sel.select(timeout=None)

        #key is the socket object 
        #mask is the event mask of operations that are ready 
        for key, mask in events: 
            accept_wrapper(key.fileobj)
        else:
            service_connection(key, mask)

except KeyboardInterrupt:
    print("Caught keyboard interrupt, exiting")

finally:
    sel.close()