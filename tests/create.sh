#!/usr/bin/env python3

from fpflow.generators.generator import Generator

generator = Generator.from_inputyaml('./input.yaml')
generator.create()
