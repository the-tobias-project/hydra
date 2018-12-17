#!/usr/bin/env python3

import _pickle as pickle


def encode(message):
    return pickle.dumps(message)


def decode(message):
    return pickle.loads(message)
