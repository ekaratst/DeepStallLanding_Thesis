import socket	
import time


s = socket.socket()		 

port = 3000				

s.bind(('', port))		 
print("socket binded to %s" %(port))

s.listen(5)	 
print("socket is listening")

c, addr = s.accept()

while True:

    c.send(str([1,2,3]).encode()) 
    time.sleep(0.2)
