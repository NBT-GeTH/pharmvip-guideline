import sys, os

def disable_print():
    sys.stdout = open(os.devnull, "w")

def enable_print():
    sys.stdout = sys.__stdout__
