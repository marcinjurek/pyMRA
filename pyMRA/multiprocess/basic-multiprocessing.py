import logging
import pdb
import multiprocessing as mp

from ipyparallel import Client



class Node(object):
    """
    The basic building block of the tree
    """
    
    
    """
    critical depth property tells you at which level in the tree you want to parallelize.
    For example if you set it to a value that's greater than the number of levels in the
    hierarchy, parallelization will not kick in (useful for debugging). The other extreme
    would be to set it to 0: then each branch of the tree starting at the root will be 
    processed separately. This might not always be desireable though: if root has, say, only
    2 children, there will be only 2 processes runnning in parallel. But if every node has 
    2 children and we set critDepth to, say, 2, then we will have 8 separate processes.
    """
    critDepth = 2
    
    def __init__(self, ID, noChildren, height):

        self.ID = ID
        self.children = []
        self.depth = len(ID)-1
        self.msg = "no message"
        if height:
            self.leaf=False
            for j in range(noChildren):
                child = Node(ID + str(j+1), noChildren, height-1)
                self.children.append( child )
        else:
            self.leaf=True

            


    def assignMsg(self, msg):
        """
        This function looks silly, but this is just a template code; In the actual
        application this function is pretty complicated
        """
        self.msg = msg

            
            
    def setMsg(self,msg, pipe=None):
        """
        This is the recursive function that sets the message for all nodes.
        Once we are at the critical depth we create as many separate processes
        as there are children, and we assign a pipe to each process. Then we send
        the child and with the "setMsg" method through the pipe for parallel
        processing. There is one caveat here though. When a child process is being
        created it becomes an exact copy of its parent (here is where the excessive
        memory usage happens). This means that when its done with processing the
        child and the subtree this child is the root of, it has to return this
        entire subtree and swap it with the existing subtree. This is not great
        in terms of efficiency but it works.
        """
        self.assignMsg(msg)
        pipes = []; procs=[]

        
        for ch in self.children:
            if self.depth==self.critDepth:
                parPipe, chPipe = mp.Pipe()
                chP = mp.Process(target=ch.setMsg, args=(msg,chPipe))
                chP.start()
                pipes.append(parPipe)
                procs.append(chP)
            else:            
                ch.setMsg(msg)


        if self.depth==self.critDepth:
            #if this is where parallelization is performed, receive the subtree
            #through the appropriate pipe and join the process.
            for idx, pipe in enumerate(pipes):
                self.children[idx] = pipe.recv()
                procs[idx].join()
                

        if self.depth==(self.critDepth+1):
            # here the subtree is returned - sent through the pipe
            pipe.send(self)

            
            

            
    def getMsg(self):
        """
        This functions returns the string that has messages from the entire
        subtree starting at this node
        """
        strOut = self.depth*'\t'+'-> %s: %s\n' % (self.ID, self.msg)
        strOut = strOut + '\n'.join([ch.getMsg() for ch in self.children])
        if self.leaf: # if it's a leaf we remove the extra newline character
            strOut = strOut[:-1]
        return strOut



    
if __name__=='__main__':
    
    MaxDepth = 2
    noChildren = 2

    root = Node('r', noChildren, MaxDepth)

    print("Before the message is set")
    print(root.getMsg())

    
    root.setMsg("I am groot")
    print("After setting the message:")
    print(root.getMsg())
