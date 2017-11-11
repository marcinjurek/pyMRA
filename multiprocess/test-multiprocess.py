import os
import traceback
import time
import numpy as np
import multiprocessing






class detCalculator(multiprocessing.Process):
  

    def __init__(self, ID, matSize, noChildren, levelsLeft, parent_pipe):

        multiprocessing.Process.__init__(self)
        self.ID = ID
        self.parent_pipe = parent_pipe
        self.mat = np.random.uniform(size=(matSize, matSize))
        self.children = []
        self.pipes = []

        print("created node " + self.ID + ":, pid=" + str(os.getpid()))

        
        if levelsLeft:
            for chNo in range(noChildren):
                chID = self.ID + str(chNo+1)
                parent_end, child_end = multiprocessing.Pipe()
                
                chP = multiprocessing.Process( target=detCalculator, args=(chID, matSize, noChildren, levelsLeft-1, child_end) )
                self.children.append( chP )
                self.pipes.append( parent_end )
                chP.start()
                




    def run(self):

        while True:
            print("Node " + self.ID + ": " + "running the loop")
            try:
                msg = self.parent_pipe.recv()
                print("Node " + self.ID + ": " + msg)
                for ch in self.children:
                    print("ch:" + ": pid=" + str(ch.pid))
                for ch_pipe in self.pipes:
                    ch_pipe.send(msg)
                    print("sending message")
            except:
                print("error")
                traceback.print_exc()
                break
        print("outside of the loop")
        return

                
            

    # def run(self):

    #     data = []
    #     while True:
    #         try:
    #             msg = self.parent_pipe.recv()
    #             if msg is None:
    #                 self.kill_children()
                    
    #             elif msg=="getDet":
    #                 data.append( (self.ID, np.linalg.det(self.mat)) )

    #                 Nremain = len(self.children)
                    
    #                 while Nremain:
    #                     if not self.q.empty():
    #                         q.get()
    #                         Nremain = Nremain - 1
    #                         data.extend(chData)

    #                 for ch_pipe in self.pipes:
    #                     ch_pipe.send(msg)

                    

                        
    #             self.parent_pipe.send( data )
    #         except EOFError:
    #             break
        







if __name__=='__main__':

    matSize = 100
    noChildren=2
    levels = 1
    
    conn1, conn2 = multiprocessing.Pipe()
    newCalc = detCalculator('r', matSize, noChildren, levels, conn2)

    newCalc.start()
    conn1.send("getDet")


    
