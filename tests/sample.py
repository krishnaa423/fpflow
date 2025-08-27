#region modules
from fpflow.inputs.grammars.siesta import SiestaGrammar
#endregion

#region variables
#endregion

#region functions
def main():
    text = r'''

SystemName       Silicon bulk
SystemLabel      Si_bulk
NumberOfAtoms    2
NumberOfSpecies  1

LatticeConstant  5.43 Ang
%block LatticeVectors
  0.5  0.5  0.0
  0.5  0.0  0.5
  0.0  0.5  0.5
%endblock LatticeVectors

%block ChemicalSpeciesLabel
  1  14  Si
%endblock ChemicalSpeciesLabel

%block AtomicCoordinatesAndAtomicSpecies
  0.0  0.0  0.0  1
  0.25 0.25 0.25 1
%endblock AtomicCoordinatesAndAtomicSpecies


MeshCutoff        200 Ry


%block kgrid_Monkhorst_Pack
  4  0  0  0
  0  4  0  0
  0  0  4  0
%endblock kgrid_Monkhorst_Pack


PAO.BasisSize     DZP
PAO.EnergyShift   0.02 Ry


SolutionMethod    diagon
MaxSCFIterations  100
DM.Tolerance      1.d-4


WriteCoorXmol     T
WriteForces       T

'''
    parser = SiestaGrammar()
    data = parser.read(text)
    output = parser.write(data)
    print(output)
#endregion

#region classes
#endregion

#region main
if __name__=='__main__':
    main()
#endregion