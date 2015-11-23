import sys
import getpass

import cutlass.iHMPSession

from . import project

class settings:
    username = 'rschwager'
    
def submit(sample_metadata):
    session = cutlass.iHMPSession(settings.username, getpass.getpass())
    

if __name__ == "__main__":
    sys.exit(submit(sys.argv[1]))
    
    
