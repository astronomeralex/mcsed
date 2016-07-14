# -*- coding: utf-8 -*-
"""
RETURNS THE LOCATION OF THE mcsed/ FOLDER.  FUNCTIONS LIKE A ENVIRONEMNTAL
VARIABLE BUT EASIER TO USE. Don't modify this.
"""
import os
def getHomeLocation():
    home = os.getcwd()
    home = home[:-3] # chop of "src"
    return home