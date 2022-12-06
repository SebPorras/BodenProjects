############################################
# Date: 2/12
# Author: Sebastian Porras   
# Aims: Create a simple echo practice server
# and client to practice socket programming  
# # Based off tutorial https://realpython.com/python-sockets/#reference       
############################################

import socket 

HOST = "127.0.1.1"
PORT = 65432

with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
    s.connect((HOST, PORT))

    s.sendall(b"Hello world")
    data = s.recv(1024)

print(f"Recieved {data!r}")