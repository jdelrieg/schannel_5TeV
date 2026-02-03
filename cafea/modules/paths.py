import os
import cafea 

pjoin = os.path.join

# This function takes as input any path (inside of cafea/cafea), and returns the absolute path
def cafea_path(path_in_repo):
    return pjoin(cafea.__path__[0], path_in_repo)


