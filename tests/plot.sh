#!/usr/bin/env python3

from fpflow.managers.manager import Manager

manager = Manager.from_inputyaml('./input.yaml')
manager.plot()
