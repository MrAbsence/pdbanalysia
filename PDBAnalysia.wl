(* ::Package:: *)

(* ::Text:: *)
(*This package will help you to analyze the structure file .pdb*)
(*Start in 2016-08.*)
(*Updates:*)
(*2020-1-2: add new option for ContactChartDistanceIntra to select the unit of distance.*)
(*2020-1-15: add asymmetric contact chart.*)
(*2020-3-23: add options in GetPDB to select the chain identifier.*)
(*2020-12-30: fix a bug in loading atom coords.*)
(*2021-02-27: add a function of rename entries in association.*)


(*TODO: make this package have a formal format*)


Print["PDBanalysia loaded."]


BeginPackage["PDBAnalysia`"]


(* ::Section:: *)
(*Usages of all public functions*)

(*This is a function to assist developers viewing source code*)
(*OpenProteinStructureWL::usage="Open this package file. The global varible \"Dropbox\" must be defined in init.m (check $UserBaseDirectory\\kernel)."*)


RenameKeyInAssociation::usage="RenameKeyInAssociation[assoc,oldKey,newKey]"


(*TODO: make usage more comprehensive*)
(*Only the functions or varibles listed here can be directly accessed in notebook*)
AAAbbrDict123::usage="An association used for converting Abbr. of Amino Acid, like G -> GLY"
AAAbbrDict321::usage="An association used for converting Abbr. of Amino Acid, like GLY -> G"

GetPDB::usage="GetPDB[directory,pdbfilename,LoadMode->\"Residues\"]
Load a PDB file, make each atom an association.
Usually, group all atoms by residues(LoadMode->\"Residues\")."
Options[GetPDB]={LoadMode->"Residues"
(*"Residues"(default): atom list grouped by residues,
"Atoms": just atom list,
"Graphic": return a wireframe Graphics3D,
"SegTree": return a hash tree: SegName\[Equal]> ResNum(Int)\[Equal]>AtomName \[Equal]> Infomation,
"ChainTree": return a hash tree: ChainName\[Equal]> ResNum(Int)\[Equal]>AtomName \[Equal]> Infomation*)};
LoadMode::usage="Residues, Atoms, Tree or Graphic.
\"Residues\"(default): atom list grouped by residues
\"Atoms\": just atom list
\"Graphic\": return a wireframe Graphics3D
\"SegTree\": return a hash tree: SegName\[Equal]> ResNum(Int)\[Equal]>AtomName \[Equal]> Infomation. You can access one atom associate by Tree[SegName][ResNum][AtomName]. This is used for getting atom index for constrains.
\"ChainTree\": return a hash tree: ChainName\[Equal]> ResNum(Int)\[Equal]>AtomName \[Equal]> Infomation*)"


PDBPickResidues::usage="PDBPickResidues[pdbResidueMode, range]
(range is the same as Take[], usually it is {m,n})
Take some residues from a pdb data to form a new one. It can be used to separate two chains in the same pdb file.
Important!! It can renumber all the residues based on the index in the new list"
Options[PDBPickResidues]={RenumberResidues->0};
RenumberResidues::usage="0 means not rewrite those nums, and any integer larger than 0 means rewritting and start from this num."


PDBResidueSeq::usage="Format output the residue sequence (with residue's ordinal) of the pdbdata."
Options[PDBResidueSeq]={ShowNumEvery->1(*Must \[GreaterEqual] 1*),RemoveDuplicateChains->True,ReturnSeqOnly->False};
ShowNumEvery::usage="Show numbers for every N residues. Must \[GreaterEqual] 1."
PDBResidues::usage="List all the types of residues in the pdbdata."
PDBAtoms::usage="All atoms are grouped by residues, return an association of these groups."
Options[PDBAtoms]={AtomType->"[CN]"(*Use regx to define which types of atoms should be selected*)};


PDBResidueContact::usage="List contacts between two residues"
Options[PDBResidueContact]={OutputDistance->False, DARRDistance->6.01};
PDBTabulateDARRContratsIntra::usage="List all possible DARR contacts in one molecule. It treat the input pdbdata as one molecule. If you want to get a full contact list in an oligomer, use full contact instead."
Options[PDBTabulateDARRContratsIntra]={IgnoreAdjacent->1, DARRNuclearRegx->"^C.*$", OutputStyle->"Contacts"(*"Graph","Distances","Contacts"*),SameTypeAAContact->False};
(*Contacts from the same type of AAs is usually useless in DARR*)
PDBTabulateDARRContratsInter::usage="This function is purely for inter-molecular interaction."
Options[PDBTabulateDARRContratsInter]={DARRNuclearRegx->"^C.*$",OutputStyle->"Contacts"(*"Graph","Distances","Contacts"*),SameTypeAAContact->False};
(**)
PDBTabulateDARRContratsFull::usage="This function is for a large oligomer system with many molecules."
Options[PDBTabulateDARRContratsFull]={IgnoreAdjacent->1, DARRNuclearRegx->"^C.*$", OutputStyle->"Contacts"(*"Graph","Distances","Contacts"*),SameTypeAAContact->True};


ContactChartSymmetric::usage="This is a fancy function to put all the contacts into a contact chart(table).\nIn the contact chart, the contact symbols are drawn at the positions in the contact list.\nThe default contact symbol is red box. You may change it to five angle star or circles.\nExample: ContactSymbol->Item[Graphics[Text[Style[\"\[FivePointedStar]\",Black,Bold,56]],ImageSize-> {50,50}],Background\[Rule]None],\nNonContactSymbol->Item[Graphics[{},Background-> None,ImageSize-> {50,50}],Background\[Rule]None]"
Options[ContactChartSymmetric]={ContactSymbol->Item[Graphics[{},Background-> Red,ImageSize-> {50,50}],Background->Red],NonContactSymbol->Item[Graphics[{},Background-> White,ImageSize-> {50,50}],Background->White]}
ContactSymbol::usage="ContactSymbols is an \"Items\" which will be drawn into the contact cells";
NonContactSymbol::usage="Non-contact symbols are drawn at all the remaining positions! Thus, blank is usually used as non-contact symbol."


ContactChartAsymmetric2Peptides::usage="This is a fancy function for put all the contacts into a contact chart(table). It allows you to create an asymmetric table from two peptides."
Options[ContactChartAsymmetric2Peptides]={ContactSymbol->Item[Graphics[{},Background-> Red,ImageSize-> {50,50}],Background->Red],NonContactSymbol->Item[Graphics[{},Background-> White,ImageSize-> {50,50}],Background->White]}


ContactChartDistanceIntra::usage="This function will generate one chart with the minimum distance between different residues"
Options[ContactChartDistanceIntra]={NucleusTypeRegx->"^C.*$",DARRDistance->6.01,DistanceUnit->"nm"}
DistanceUnit::usage="It only support \"nm\" or \"A\"."


ContactChartDistanceInter::usage="This function will generate one chart with the minimum distance between residues from 2 peptide molecules"
Options[ContactChartDistanceInter]={NucleusTypeRegx->"^C.*$",DARRDistance->6.01,DistanceUnit->"nm"}


Begin["`Private`"](*The following is private functions*)


(* ::Section:: *)
(*Some common knowledge of this package*)


(*These two associations are used for converting Abbr. of Amino Acid*)
AAAbbrDict123=<|"F"->"PHE","L"->"LEU","S"->"SER","Y"->"TYR","C"->"CYS","W"->"TRP","P"->"PRO","H"->"HIS","Q"->"GLN","R"->"ARG","I"->"ILE","M"->"MET","T"->"THR","N"->"ASN","K"->"LYS","V"->"VAL","A"->"ALA","D"->"ASP","E"->"GLU","G"->"GLY"|>;
AAAbbrDict321=Association[Reverse/@Normal[AAAbbrDict123]];

(*A private association referencing the entries in an atom association*)
(*An atom must have an association like this. It's a reference for all the entries*)
AtomAssociasionExample=<|"AtomNums"->"Int","ElementType"->"1Letter","AtomNames"->"Str","ResidueNums"->"Int","ResidueNames"->"3Letter","ChainIDs"->"1Letter","AtomCoords"->"3Floats","SegmentName"->"Str"|>;


(*Some private functions come from some dependences.*)
(*ChangeDirctory is from NMRSpectrumAnalysis*)
ChangeDirectory[dirstring_]:=If[$OperatingSystem=="Unix" || $OperatingSystem == "MacOSX",SetDirectory[StringReplace[dirstring,"\\"->"/"]],SetDirectory[StringReplace[dirstring,"/"->"\\\\"]]]


RenameKeyInAssociation[assoc_,oldKey_,newKey_]:=Module[{assocNew},
assocNew = assoc;
assocNew[newKey] = assoc[oldKey];
assocNew = KeyDrop[assocNew,oldKey];
assocNew]


(*supporting functions that might be useful some day*)

(*A function to open the source code of this package*)
(*OpenProteinStructureWL[]:=Button["open ProteinStructure.wl",SystemOpen[Global`Dropbox<>"\\MathematicaConfig\\YG\\ProteinStructure.wl"]]*)


(* ::Section:: *)
(*Load the PDB information*)


(*Load a PDB file, group all atoms by residues*)
GetPDB[wd_, pdbfile_,opts___?OptionQ] :=Module[{PDBData,PDBName,AtomCoords,AtomNums,AtomNames,ResidueNames,ResidueNums,ElementType,ChainIDs,SegmentName,lm},
lm=LoadMode/. {opts}/. Options[GetPDB];
ChangeDirectory[wd];
(*append the file name with .pdb*)
PDBName=If[StringFreeQ[pdbfile,".pdb"],pdbfile<>".pdb",pdbfile];
(*Some options will determine how the file is loading*)
Which[lm=="Graphic",
(*Graphic: get a wireframe Graphics3D*)
Print["Load as 3D graph"];
PDBData=Import[PDBName,"Rendering"->"Wireframe"],
(*A trick for Which, True at last fits all other conditions*)
True,
(*Grab all atom lines into list*)
(*Import the PDB file by lines*)
PDBData=Import[PDBName,"Lines"];
(*Select all atom lines from PDB file*)
PDBData =Select[PDBData,StringMatchQ[#,RegularExpression["^ATOM\\s+\\d+.+$"]]&];
(*Get useful informations in to different lists*)
(*Based on PDB format: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/framepdbintro.html*)
(*The string split method may work in most conditions. However, there is no space between those numbers in some pdb file. e.g. 2LNQ*)
(*AtomCoords =ToExpression[#]&/@StringSplit[StringTake[PDBData,{31,54}]];*)
AtomCoords = ToExpression[#]&/@StringTake[PDBData,{{31,38},{39,46},{47,54}}];
AtomNums =ToExpression[#]&/@StringTake[PDBData,{6,11}];
ResidueNums=ToExpression[#]&/@StringTake[PDBData,{23,26}];
AtomNames = StringTrim[#]&/@StringTake[PDBData,{13,16}];
ElementType=StringTake[StringTrim[#],1]&/@StringTake[PDBData,{13,14}];
ResidueNames=StringTake[PDBData,{18,20}];
ChainIDs=StringTake[PDBData,{22}];
SegmentName=StringTrim[#]&/@StringTake[PDBData,{73,76}];
(*Segment identifier is obsolete, but still used by some programs. Chimera assigns it as the atom attribute pdbSegment to allow command-line specification.*)
(*Format the atoms' info into one association per atom. Use AtomAssociasionExample as ref for entries*)
PDBData=Transpose[{AtomNums,ElementType,AtomNames,ResidueNums,ResidueNames,ChainIDs,AtomCoords,SegmentName}];
PDBData=AssociationThread[Keys[AtomAssociasionExample]->#]&/@PDBData;
(*Split the atom list by residues or not split*)
PDBData=Which[lm=="Residues",
(*If there are two adjecent residues named G280 and A280, this would go wrong*)
SplitBy[PDBData,#ResidueNums&],
lm=="SegTree",
(*SegName-ResID-AtomName Tree*)
Map[GroupBy[#,#["AtomNames"]&,First]&,(GroupBy[#,#["ResidueNums"]&]&)/@GroupBy[PDBData,#["SegmentName"]&],{2}],
lm=="ChainTree",
Map[GroupBy[#,#["AtomNames"]&,First]&,(GroupBy[#,#["ResidueNums"]&]&)/@GroupBy[PDBData,#["ChainIDs"]&],{2}],
lm=="Atoms",
PDBData
];
(*Print the number of residues or atoms for reference*)
Which[lm=="Residues",
Print[ToString[Length[PDBData]]<>" Residues"],
lm=="Atoms",
Print[ToString[Length[PDBData]]<>" Atoms"]]
];
ResetDirectory[];
PDBData
]

(*Take some residues from a pdb data to form a new one. It could be used to separate two chains in the same pdb file. Important!! It can renumber all the residues based on the index in the new list*)
PDBPickResidues[pdbdata_,residualnum_,opts___?OptionQ]:=
Module[{PickedResidues,rnr},
rnr=RenumberResidues/.{opts}/.Options[PDBPickResidues];
PickedResidues=Take[pdbdata,residualnum];
If[rnr==0,PickedResidues,Table[<|#,"ResidueNums"->i+rnr-1|>&/@PickedResidues[[i]],{i,1,Length[PickedResidues]}]]
]


(* ::Section:: *)
(*Manipulate content in PDB files*)


PDBFileAssignChains[pdbfile_,numberOfAtomsEachChain_,opts___?OptionQ]:=Module[{pdblines,
chainNameList=Characters["ABCDEFGHIJKLMNOPQRSTUVWXYZ"]},
pdblines=Import[pdbfile,"Line"];
pdblines(*not done*)
]


(* ::Section:: *)
(*Display residue info in PDBData*)


(*Format output the residue sequence (with residue's ordinal) of the pdbdata*)
PDBResidueSeq[pdbdata_,opts___?OptionQ]:=
Module[{AASeq,sne,rmdu,rso},
sne=ShowNumEvery/. {opts}/. Options[PDBResidueSeq];
rmdu=RemoveDuplicateChains/.{opts}/.Options[PDBResidueSeq];
rso=ReturnSeqOnly/.{opts}/.Options[PDBResidueSeq];
AASeq=AAAbbrDict321[#[[1]]["ResidueNames"]]<>ToString[#[[1]]["ResidueNums"]]&/@pdbdata;
(*If remove duplicated chains, just remove identical elements in AASeq; If not remove duplicated chains, label all the duplicated AAs*)
AASeq=If[rmdu,DeleteDuplicates[AASeq],Table[If[Count[AASeq[[1;;i]],AASeq[[i]]]>1,AASeq[[i]]<>"-"<>ToString[Count[AASeq[[1;;i]],AASeq[[i]]]],AASeq[[i]]],{i,Length[AASeq]}]];
(*This can only work on duplicated chains, not different chains*)
If[rso,AASeq,
Style[StringJoin[Table[If[Mod[i-1,sne]==0,AASeq[[i]]<>",",StringTake[AASeq[[i]],1]],{i,Length[AASeq]}]],Blue,Bold,20]]
]
(*ReturnSeqOnly: return a list of all the residues. If false, then return a format string.*)


(*List all the types of residues in the pdbdata*)
PDBResidues[pdbdata_]:=
Module[{},Union[AAAbbrDict321[#]&/@(#[[1]]["ResidueNames"]&/@pdbdata)]]


(*All atoms are grouped by residues, return an association of these groups*)
PDBAtoms[pdbdata_,opts___?OptionQ]:=
Module[{at,PickedAtoms},
at=AtomType/.{opts}/.Options[PDBAtoms];
(*Pick atoms based on ElementType*)
PickedAtoms=Select[#,StringMatchQ[#ElementType,RegularExpression[atomType]]&]&/@pdbdata;
AssociationThread[PDBResidueSeq[pdbdata,RemoveDuplicateChains->False,ReturnSeqOnly->True],Map[#AtomNames&,PickedAtoms,{2}]]
]



(* ::Section:: *)
(*Detect contacts between residues*)


(*Inputs must be lists consist of the same type of atoms*)
PDBResidueContact[residueatoms1_,residueatoms2_,opts___?OptionQ]:=
Module[{od,ddis,DistanceTable,FormatDistanceTable,FormatAssist},
od=OutputDistance/.{opts}/.Options[PDBResidueContact];
ddis=DARRDistance/.{opts}/.Options[PDBResidueContact];
(*Calculate the distances between each atoms*)
DistanceTable=Table[EuclideanDistance[residueatoms1[[i]]["AtomCoords"],residueatoms2[[j]]["AtomCoords"]],{i,1,Length[residueatoms1]},{j,1,Length[residueatoms2]}];
(*Format the output table*)
FormatDistanceTable=Insert[DistanceTable,AAAbbrDict321[#ResidueNames]<>ToString[#ResidueNums]<>#AtomNames&/@residueatoms2,1];
FormatAssist=Insert[AAAbbrDict321[#ResidueNames]<>ToString[#ResidueNums]<>#AtomNames&/@residueatoms1,AAAbbrDict321[residueatoms1[[1]]["ResidueNames"]]<>ToString[residueatoms1[[1]]["ResidueNums"]]<>"\\"<>AAAbbrDict321[residueatoms2[[1]]["ResidueNames"]]<>ToString[residueatoms2[[1]]["ResidueNums"]],1];
FormatDistanceTable=Insert[Transpose[FormatDistanceTable],FormatAssist,1];
(*If not printing distance, replace distance larger than DARRDistance with "No Cnt"*)
If[od,FormatDistanceTable,FormatDistanceTable/.dis_?NumericQ/;dis>ddis->"No Cnt"] (*Distance>6A \[Equal]> No Contact!*)
]


(*Only return the shortest distance between the two residue (same type of nucleus, usually 13C)*)
PDBResidueMinDistance[residueatoms1_,residueatoms2_,opts___?OptionQ]:=
Module[{DistanceMin},
(*Calculate the distances between each atoms*)
DistanceMin=Min[Table[EuclideanDistance[residueatoms1[[i]]["AtomCoords"],residueatoms2[[j]]["AtomCoords"]],{i,1,Length[residueatoms1]},{j,1,Length[residueatoms2]}]]
]


PDBTabulateDARRContratsIntra[pdbdata_,opts___?OptionQ]:=
Module[{igad,dnr,osty,stac,AAs,DistanceGroups,ContactTables},
igad=IgnoreAdjacent/.{opts}/.Options[PDBTabulateDARRContratsIntra];
dnr=DARRNuclearRegx/.{opts}/.Options[PDBTabulateDARRContratsIntra];(*In fact this is not necessary, we can use the AtomNames*)
osty=OutputStyle/.{opts}/.Options[PDBTabulateDARRContratsIntra];
stac=SameTypeAAContact/.{opts}/.Options[PDBTabulateDARRContratsIntra];
(*Selecting only the same type of nucleus by DARRNuclearRegx*)
(*The ResidueNums\[Rule]ResAtomInfo format will help to store the residue sequence so that we can make sure what is adjacent*)
AAs=Association@@((#[[1]]["ResidueNums"]->Select[#,StringMatchQ[#AtomNames,RegularExpression[dnr]]&])&/@pdbdata);
(*List all possible groups need to test*)
DistanceGroups=Select[Subsets[Keys[AAs],{2}],Abs[Subtract[#[[1]],#[[2]]]]>igad&];
(*Calculate contact or not*)
ContactTables=PDBResidueContact[AAs[#[[1]]],AAs[#[[2]]]]&/@DistanceGroups;
(*Refine for contacts: based on the return from PDBResidueContact, any number \[Equal]> contact*)
ContactTables=Select[ContactTables,AnyTrue[#,NumericQ,2]&];
(*Eliminate Same type AA contacts*)
ContactTables=If[stac,ContactTables,Select[ContactTables,!StringMatchQ[Part[#,1,1],RegularExpression["^([A-Z])\\d+\\\\\\1\\d+$"]]&]];
(*Different types of ouput*)
Which[osty=="Graph",
GraphPlot[ToExpression[StringReplace[Part[#,1,1],"\\"->"->"]&/@ContactTables],VertexLabeling->True],osty=="Contacts",Part[#,1,1]&/@ContactTables,
osty=="Distances",ContactTables//TableForm]
]



(*This function is purely for inter-molecular interactions*)
PDBTabulateDARRContratsInter[pdbdata1_,pdbdata2_,opts___?OptionQ]:=
Module[{dnr,osty,stac,DistanceGroups,ContactTables,AAs1,AAs2},
dnr=DARRNuclearRegx/.{opts}/.Options[PDBTabulateDARRContratsInter];(*In fact this is not necessary, we can use the AtomNames*)
osty=OutputStyle/.{opts}/.Options[PDBTabulateDARRContratsInter];
stac=SameTypeAAContact/.{opts}/.Options[PDBTabulateDARRContratsInter];
Print["Chain 1:"<>PDBResidueSeq[pdbdata1][[1]]];
Print["Chain 2:"<>PDBResidueSeq[pdbdata2][[1]]];
(*Select specific nuclear by regx*)
AAs1=Select[#,StringMatchQ[#AtomNames,RegularExpression[dnr]]&]&/@pdbdata1;
AAs2=Select[#,StringMatchQ[#AtomNames,RegularExpression[dnr]]&]&/@pdbdata2;
ContactTables=Flatten[Table[PDBResidueContact[AAs1[[i]],AAs2[[j]]],{i,Length[AAs1]},{j,Length[AAs2]}],1];
(*Refine for contacts: based on the return from PDBResidueContact, any number \[Equal]> contact*)
ContactTables=Select[ContactTables,AnyTrue[#,NumericQ,2]&];
(*Eliminate Same type AA contacts*)
ContactTables=If[stac,ContactTables,Select[ContactTables,!StringMatchQ[Part[#,1,1],RegularExpression["^([A-Z])\\d+\\\\\\1\\d+$"]]&]];
(*Different types of ouput*)
Which[osty=="Graph",
GraphPlot[ToExpression[StringReplace[Part[#,1,1]<>"-2","\\"->"->"]&/@ContactTables],VertexLabeling->True],osty=="Contacts",(Part[#,1,1]<>"-2")&/@ContactTables,
osty=="Distances",ContactTables//TableForm]
]


PDBTabulateDARRContratsFull[pdbdata_,opts___?OptionQ]:=
Module[{igad,dnr,osty,stac,AAs,DistanceGroups,ContactTables},
igad=IgnoreAdjacent/.{opts}/.Options[PDBTabulateDARRContratsFull];
dnr=DARRNuclearRegx/.{opts}/.Options[PDBTabulateDARRContratsFull];(*In fact this is not necessary, we can use the AtomNames*)
osty=OutputStyle/.{opts}/.Options[PDBTabulateDARRContratsFull];(*For extend in the future.*)
stac=SameTypeAAContact/.{opts}/.Options[PDBTabulateDARRContratsFull];
(*Selecting only the same type of nucleus by DARRNuclearRegx*)
(*The ResidueNums\[Rule]ResAtomInfo format will help to store the residue sequence so that we can make sure what is adjacent*)
AAs=Select[#,StringMatchQ[#AtomNames,RegularExpression[dnr]]&]&/@pdbdata;
(*List all possible groups need to test*)
DistanceGroups=Select[Subsets[AAs,{2}],Abs[Subtract[#[[1,1]]["ResidueNums"],#[[2,1]]["ResidueNums"]]]>igad&];
(*Calculate contact or not*)
ContactTables=PDBResidueContact[#[[1]],#[[2]]]&/@DistanceGroups;
(*Refine for contacts: based on the return from PDBResidueContact, any number \[Equal]> contact*)
ContactTables=Select[ContactTables,AnyTrue[#,NumericQ,2]&];
(*Eliminate Same type AA contacts*)
ContactTables=If[stac,ContactTables,Select[ContactTables,!StringMatchQ[Part[#,1,1],RegularExpression["^([A-Z])\\d+\\\\\\1\\d+$"]]&]];
(*Different types of ouput*)
Which[osty=="Contacts",Part[#,1,1]&/@ContactTables,
osty=="FUTUREFUNCTION",0]
]


(* ::Section:: *)
(*Display Functions*)


(* ::Input:: *)
(*(*TODO: Add functions for making figures!*)*)


(*This is a fancy function for put all the contacts into a contact chart(table)*)
(*In the contact chart, the contact symbols are drawn at the positions in the contact list.*)
(*Non-contact symbols are drawn at all the remaining positions! Thus, blank is usually used as non-contact symbol.*)
(*For example, if "A1\F10" is in the contact list. There will be a contact symbol at the cross cell of A1 and F10.*)
ContactChartSymmetric[aalist_,contactlist_,opts___?OptionQ]:=
Module[{csys,ncsys,indexcontactlist,contactsyltable,fulltable},
csys=ContactSymbol/.{opts}/.Options[ContactChartSymmetric];
ncsys=NonContactSymbol/.{opts}/.Options[ContactChartSymmetric];
indexcontactlist=Flatten@(Position[aalist,#]&/@(StringSplit[#,"\\"]))&/@contactlist;
indexcontactlist=Flatten[{#,Reverse[#]}&/@indexcontactlist,1];
indexcontactlist=Union[indexcontactlist];
contactsyltable=Table[If[MemberQ[indexcontactlist,{ResIdx1,ResIdx2}],csys,ncsys],{ResIdx1,1,Length[aalist]},{ResIdx2,1,Length[aalist]}];
fulltable=Prepend[contactsyltable,Item[Graphics[Text[Style[#,Bold,28,FontFamily->"Arial"]],ImageSize->{50,50}],Background->None]&/@aalist];
fulltable=Insert[Transpose [fulltable],Prepend[Item[Graphics[Text[Style[#,Bold,28,FontFamily->"Arial"]],ImageSize->{50,50}],Background->None]&/@aalist, ""],1]//Transpose; (*Add extra grid box to aalist so matrix dimensions match*)
Grid[fulltable, Frame-> All,Spacings->{1.5,1.5},ItemSize->Full]
]


ContactChartAsymmetric2Peptides[aalist1_,aalist2_,contactlist_,opts___?OptionQ]:=
Module[{csys,ncsys,indexcontactlist,contactsyltable,fulltable},
csys=ContactSymbol/.{opts}/.Options[ContactChartSymmetric];
ncsys=NonContactSymbol/.{opts}/.Options[ContactChartSymmetric];
indexcontactlist={First@@Position[aalist1,#[[1]]],First@@Position[aalist2,#[[2]]]}&/@((StringSplit[#,"\\"])&/@contactlist);
contactsyltable=Table[If[MemberQ[indexcontactlist,{ResIdx1,ResIdx2}],csys,ncsys],{ResIdx1,1,Length[aalist1]},{ResIdx2,1,Length[aalist2]}];
fulltable=Prepend[contactsyltable,Item[Graphics[Text[Style[#,Bold,28,FontFamily->"Arial"]],ImageSize->{50,50}],Background->None]&/@aalist2];
fulltable=Insert[Transpose [fulltable],Prepend[Item[Graphics[Text[Style[#,Bold,28,FontFamily->"Arial"]],ImageSize->{50,50}],Background->None]&/@aalist1, ""],1]//Transpose; (*Add extra grid box to aalist so matrix dimensions match*)
Grid[fulltable, Frame-> All,Spacings->{1.5,1.5},ItemSize->Full]
]


(*This function will generate one chart with the minimum distance between different residues in one peptide*)
ContactChartDistanceIntra[pdbdata_,opts___?OptionQ]:=
Module[{ntype,ddis,aalist, filterpdbdata,distanceTable,fulltable,di,dunit},
ntype=NucleusTypeRegx/.{opts}/.Options[ContactChartDistanceIntra];
ddis=DARRDistance/.{opts}/.Options[ContactChartDistanceIntra];
dunit=DistanceUnit/.{opts}/.Options[ContactChartDistanceIntra];
aalist=PDBResidueSeq[pdbdata,ReturnSeqOnly->True];
(*Filter the nucleus we care in the pdbdata*)
filterpdbdata=(Select[#,StringMatchQ[#AtomNames,RegularExpression[ntype]]&]&)/@pdbdata;
distanceTable=Table[Item[Graphics[Text[Style[
ToString[NumberForm[If[dunit=="A",di=PDBResidueMinDistance[filterpdbdata[[ResIdx1]],filterpdbdata[[ResIdx2]]](*A*),(di=PDBResidueMinDistance[filterpdbdata[[ResIdx1]],filterpdbdata[[ResIdx2]]])/10.0(*nm*)],{3,2}]]
,Bold,20,If[di<ddis,White,Black],FontFamily->"Arial"]],ImageSize->{50,50}],Background->If[di<ddis,Gray,None]]
,{ResIdx1,1,Length[aalist]},{ResIdx2,1,Length[aalist]}];
fulltable=Prepend[distanceTable,Item[Graphics[Text[Style[#,Bold,28,FontFamily->"Arial"]],ImageSize->{50,50}],Background->None]&/@aalist];
fulltable=Insert[Transpose [fulltable],Prepend[Item[Graphics[Text[Style[#,Bold,28,FontFamily->"Arial"]],ImageSize->{50,50}],Background->None]&/@aalist, ""],1]//Transpose; (*Add extra grid box to aalist so matrix dimensions match*)
Grid[fulltable, Frame-> All,Spacings->{1.5,1.5},ItemSize->Full]
]


(*This function will generate one chart with the minimum distance between residues from 2 peptide molecules*)
ContactChartDistanceInter[pdbdata1_,pdbdata2_,opts___?OptionQ]:=
Module[{ntype,ddis,aalist1,aalist2,filterpdbdata1,filterpdbdata2,distanceTable,fulltable,di,dunit},
ntype=NucleusTypeRegx/.{opts}/.Options[ContactChartDistanceInter];
ddis=DARRDistance/.{opts}/.Options[ContactChartDistanceInter];
dunit=DistanceUnit/.{opts}/.Options[ContactChartDistanceIntra];
aalist1=PDBResidueSeq[pdbdata1,ReturnSeqOnly->True];
aalist2=PDBResidueSeq[pdbdata2,ReturnSeqOnly->True];
(*Filter the nucleus we care in the pdbdata*)
filterpdbdata1=(Select[#,StringMatchQ[#AtomNames,RegularExpression[ntype]]&]&)/@pdbdata1;
filterpdbdata2=(Select[#,StringMatchQ[#AtomNames,RegularExpression[ntype]]&]&)/@pdbdata2;
distanceTable=Table[Item[Graphics[Text[Style[
ToString[NumberForm[If[dunit=="A",di=PDBResidueMinDistance[filterpdbdata[[ResIdx1]],filterpdbdata[[ResIdx2]]](*A*),(di=PDBResidueMinDistance[filterpdbdata[[ResIdx1]],filterpdbdata[[ResIdx2]]])/10.0(*nm*)],{4,2}]]
,Bold,20,If[di<ddis,White,Black],FontFamily->"Arial"]],ImageSize->{50,50}],Background->If[di<ddis,Red,None]]
,{ResIdx1,1,Length[aalist1]},{ResIdx2,1,Length[aalist2]}];
fulltable=Prepend[distanceTable,Item[Graphics[Text[Style[#,Bold,28,FontFamily->"Arial"]],ImageSize->{50,50}],Background->None]&/@aalist2];
fulltable=Insert[Transpose [fulltable],Prepend[Item[Graphics[Text[Style[#,Bold,28,FontFamily->"Arial"]],ImageSize->{50,50}],Background->None]&/@aalist1, ""],1]//Transpose; (*Add extra grid box to aalist so matrix dimensions match*)
Grid[fulltable, Frame-> All,Spacings->{1.5,1.5},ItemSize->Full]
]


End[]


EndPackage[]
