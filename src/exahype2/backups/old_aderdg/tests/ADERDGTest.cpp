// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#include "ADERDGTest.h"

#include <sstream>
#include <iostream>
#include <iomanip>

#include <limits>
#include <cmath>

#include "tarch/la/Vector.h"

#include "exahype2/backups/old_aderdg/KernelUtils.h"
#include "exahype2/backups/old_aderdg/PredictorAoS.h"
#include "exahype2/backups/old_aderdg/RiemannAoS.h"
#include "exahype2/backups/old_aderdg/RusanovNonlinearAoS.h"
#include "exahype2/backups/old_aderdg/CorrectorAoS.h"

exahype2::aderdg::tests::ADERDGTest::ADERDGTest():
  TestCase ("exahype2::aderdg::ADERDGTest"),
  QuadraturePoints{0.069431844202973713731097404888714663684368133544921875,0.33000947820757187134432797392946667969226837158203125,0.66999052179242812865567202607053332030773162841796875,0.9305681557970262307577513638534583151340484619140625},
  QuadratureWeights{0.1739274225687269248563637802362791262567043304443359375,0.32607257743127304738806060413480736315250396728515625,0.32607257743127304738806060413480736315250396728515625,0.1739274225687269248563637802362791262567043304443359375},
  BarycentricWeights{-7.4205400680389477230392003548331558704376220703125,18.79544940755506132745722425170242786407470703125,-18.79544940755506132745722425170242786407470703125,7.42054006803894861121762005495838820934295654296875},
  BasisFunctionValuesLeft{1.5267881254572663873858573424513451755046844482421875,-0.8136324494869271450880887641687877476215362548828125,0.40076152031165046540905905203544534742832183837890625,-0.11391719628198997138479597879268112592399120330810546875},
  BasisFunctionValuesRight{-0.113917196281990040773735017864964902400970458984375,0.400761520311650742964815208324580453336238861083984375,-0.813632449486927811221903539262712001800537109375,1.5267881254572674976088819676078855991363525390625},
  DerivativeOperator{-6.6640004727045631938153746887110173702239990234375,9.7203088313703940315235740854404866695404052734375,-4.21756469699035818621268845163285732269287109375,1.161256338324529568950538305216468870639801025390625,-1.5151152295984677831341969067580066621303558349609375,-0.7688287844464174458636307463166303932666778564453125,2.941340462561433444221847821609117090702056884765625,-0.6573964485165488813578349436284042894840240478515625,0.65739644851654832624632263105013407766819000244140625,-2.941340462561433444221847821609117090702056884765625,0.768828784446416335640606121160089969635009765625,1.515115229598468449268011681851930916309356689453125,-1.161256338324528680772118605091236531734466552734375,4.21756469699035729803426875150762498378753662109375,-9.720308831370392255166734685190021991729736328125,6.66400047270456497017221408896148204803466796875},
  MassMatrix{0.1739274225687269248563637802362791262567043304443359375,0.0,0.0,0.0,0.0,0.32607257743127304738806060413480736315250396728515625,0.0,0.0,0.0,0.0,0.32607257743127304738806060413480736315250396728515625,0.0,0.0,0.0,0.0,0.1739274225687269248563637802362791262567043304443359375},
  StiffnessOperator{-1.1590524262142825051569161587394773960113525390625,-0.494037528020547400675610560938366688787937164306640625,0.214358954361956122180998818294028751552104949951171875,-0.201974321866382811041518152705975808203220367431640625,1.69062826161228674237690938753075897693634033203125,-0.250693983347795967819848783619818277657032012939453125,-0.95909046573029943516530693159438669681549072265625,0.73355015726438665968345276269246824085712432861328125,-0.7335501572643867707057552252081222832202911376953125,0.95909046573029943516530693159438669681549072265625,0.250693983347795634752941396072856150567531585693359375,-1.690628261612286298287699537468142807483673095703125,0.2019743218663829775749718464794568717479705810546875,-0.21435895436195628871445251206750981509685516357421875,0.494037528020547622720215485969674773514270782470703125,1.159052426214282949246126008802093565464019775390625},
  InvertedPredictorLhsOperator{0.546435362419644660267499525208023669620160944759845733642578125,-0.14432618329325673533526652736469486626447178423404693603515625,0.1014623598288631781773111092959105405952868750318884849548339844,-0.06687578310368458280271795196592066190532932523638010025024414063,1.0188533167713026317442726043083212061901576817035675048828125,0.5847599728573232642954426996340089317527599632740020751953125,-0.1702637242678432449784643037959952494020399171859025955200195313,0.1014623598288624077161422867843221240491402568295598030090332031,1.0240105066930918847022125017787175238481722772121429443359375,1.0007437785531960714770216558378024274134077131748199462890625,0.58475997285732400610656911421614267965196631848812103271484375,-0.1443261832932565697233846802038925716260564513504505157470703125,0.97400505826439649705202061813480440832790918648242950439453125,1.02401050669309141394162920857269227781216613948345184326171875,1.0188533167713030021077347253566358631360344588756561279296875,0.546435362419644295325048266587231182711548171937465667724609375},
  EquidistantGridProjector{1.5267881254572663873858573424513451755046844482421875,-0.0049591884824765342099084364235750399529933929443359375,-0.0021913277260105974188209021491502426215447485446929931640625,-0.113917196281990040773735017864964902400970458984375,-0.8136324494869271450880887641687877476215362548828125,0.997304018973186767738070557243190705776214599609375,0.0098464972353006600946923043693459476344287395477294921875,0.400761520311650742964815208324580453336238861083984375,0.40076152031165046540905905203544534742832183837890625,0.00984649723530049529596208657267197850160300731658935546875,0.99730401897318643467116316969622857868671417236328125,-0.813632449486927811221903539262712001800537109375,-0.11391719628198997138479597879268112592399120330810546875,-0.00219132772601056229067051361880658078007400035858154296875,-0.004959188482476617476635283310315571725368499755859375,1.5267881254572674976088819676078855991363525390625},
  FineGridProjector{{1.3365809435384148340375531915924511849880218505859375,0.75017378205554174908087361473008058965206146240234375,0.25006831351327185597455127208377234637737274169921875,0.03282922721362747930928804862560355104506015777587890625,-0.5106599509271045889136075857095420360565185546875,0.350399125663716670686653742450289428234100341796875,0.9137546431442424843538674394949339330196380615234375,1.01007139472141194147525311564095318317413330078125,0.2422582771987365768406874622087343595921993255615234375,-0.1376638598070825392216676164025557227432727813720703125,-0.2182390043731604889476471953457803465425968170166015625,-0.055641038587549929150810612554778344929218292236328125,-0.068179269810046794209057452462729997932910919189453125,0.037090952087824181904185394387241103686392307281494140625,0.054416047715646100046971156416475423611700534820556640625,0.01274041665251061418440148287345436983741819858551025390625},{-0.035350042596419377349814539002181845717132091522216796875,-0.0928683884767426970352488524440559558570384979248046875,-0.071267786528997401074292383782449178397655487060546875,-0.01767502129820965051099079801133484579622745513916015625,0.97104619867609065497759956997469998896121978759765625,0.7760907833371601949323803637525998055934906005859375,0.388045391668580041955038950618472881615161895751953125,0.08197886521853813002191913028582348488271236419677734375,0.081978865218538310433160631873761303722858428955078125,0.38804539166857987542158525684499181807041168212890625,0.776090783337159972887775438721291720867156982421875,0.97104619867609087702220449500600807368755340576171875,-0.0176750212982097025526950773155476781539618968963623046875,-0.07126778652899741495208019159690593369305133819580078125,-0.0928683884767427525464000837018829770386219024658203125,-0.035350042596419335716451115558811579830944538116455078125},{0.01274041665251053785656853989394221571274101734161376953125,0.054416047715646016780244309529734891839325428009033203125,0.037090952087824084759670739686043816618621349334716796875,-0.0681792698100467664534818368338164873421192169189453125,-0.055641038587549630778372744543958106078207492828369140625,-0.2182390043731603779253447328301263041794300079345703125,-0.1376638598070822616659114601134206168353557586669921875,0.242258277198736549085111846579820849001407623291015625,1.01007139472141194147525311564095318317413330078125,0.9137546431442427063984723645262420177459716796875,0.350399125663715838019385273582884110510349273681640625,-0.5106599509271045889136075857095420360565185546875,0.032829227213627291959152643130437354557216167449951171875,0.250068313513271689441097578310291282832622528076171875,0.7501737820555423041923859273083508014678955078125,1.3365809435384148340375531915924511849880218505859375}},
  testEuler_UIn{
    1.00000000000000,  0.100000000000000, 0.200000000000000, 0.300000000000000, 2.57000000000000},
  testEuler_QOut{
    1.00000000000000,  0.100000000000000, 0.200000000000000, 0.300000000000000, 2.57000000000000},
  testAdvection_UIn{
    1.00000000000000,  1.000000000000000, 1.123456789000000, 0.987654321000000}
{
}

void exahype2::aderdg::tests::ADERDGTest::runADERDGStep(
  std::function< void(
    const double * const __restrict__           Q,
    const tarch::la::Vector<Dimensions,double>& x,
    double                                      t,
    int                                         normal,
    double * __restrict__                       F
  ) >                                          flux,
  std::function< double(
    const double * const __restrict__           Q,
    const tarch::la::Vector<Dimensions,double>& x,
    double                                      t,
    int                                         normal
  ) >                                          maxEigenvalue,
  std::function< void(
    const double* Q,
    const int     i_Q,
    const int     unknown
  )>                                           validatePicardLoopHullAndCorrectorResult,
  const tarch::la::Vector<Dimensions, double>& x,
  const double                                 dx,
  const double                                 t,
  const double                                 dt,
  const int                                    unknowns,
  const int                                    auxiliaryVariables,
  const int                                    order,
  const double*                                test_UIn,
  const bool                                   verbose
) {
  std::ostringstream buffer;

  const int strideQ       = unknowns+auxiliaryVariables;
  const int nodesPerAxis  = (order+1);
  int nodesPerFace  = nodesPerAxis; 
  if ( Dimensions == 3 ) {
    nodesPerFace *= nodesPerAxis;
  }
  const int nodesPerCell          = nodesPerFace*nodesPerAxis;
  const int spaceTimeNodesPerCell = nodesPerCell*nodesPerAxis;

  // In-/Outputs:
  buffer << "\ninput (U):\n\n";
  double U[nodesPerCell*strideQ];
  for (int i = 0; i < nodesPerCell; i++) {
    for (int m=0; m < unknowns; m++) {
      const int i_U = i*strideQ + m;
      U[i_U] = test_UIn[m];
      buffer << std::fixed << std::showpoint << std::setprecision(6) << U[i_U] << " ";
    }
    for (int m=unknowns; m < strideQ; m++) {
      const int i_U = i*strideQ + m;
      U[i_U] = 0.0;
      buffer << std::fixed << std::showpoint << std::setprecision(6) << U[i_U] << " ";
    }
    buffer << "\n";
  }
  if ( verbose ) {
    std::cout << buffer.str() << std::flush; 
  }

  // Outputs:
  double Q[spaceTimeNodesPerCell*strideQ];
  double _QHull[Dimensions*2][2*nodesPerCell*strideQ];
  double* QHull[Dimensions*2];
  double* QHull_faceSwap[Dimensions*2]; // left becomes right
  double riemannResult[Dimensions*2*nodesPerAxis];
  bool atBoundary[Dimensions*2];
  for (int face = 0; face < Dimensions*2; face++) {
    const int orientation = face/Dimensions;
    const int direction   = face - Dimensions*orientation;
    QHull         [ face ] = &_QHull[face][0];
    QHull_faceSwap[ face ] = &_QHull[Dimensions*(1-orientation)+direction][0];
    atBoundary    [ face ] = false;
  }
  double maxEigenvaluePerFace[Dimensions*2];

  ::exahype2::aderdg::spaceTimePredictor_PicardLoop_loop_AoS(
    flux,
    [&](
      const double * const __restrict__           Q,
      const tarch::la::Vector<Dimensions,double>& x,
      double                                      t,
      double * __restrict__                       S
    )->void {},
    [&](
      const double * const __restrict__           Q,
      const double * const __restrict__           dQ_or_dQdn,
      const tarch::la::Vector<Dimensions,double>& x,
      double                                      t,
      int                                         direction,
      double * __restrict__                       BgradQ
    )->void {},
    Q,                                 // Q
    U,                                 // U
    QuadraturePoints,
    QuadratureWeights,
    StiffnessOperator,                 // Kxi,
    InvertedPredictorLhsOperator,      // iK1,
    BasisFunctionValuesLeft,           // FLCoeff,
    DerivativeOperator,                // dudx, 
    x,
    dx,                                // we assume cubic/square cells
    t,
    dt,
    order, 
    unknowns,
    auxiliaryVariables,
    1e-8,  // atol,
    true, 
    false,
    false
  );
  
  // print result
  buffer << "\nresult (Q):\n\n";
  for (int i = 0; i < spaceTimeNodesPerCell; i++) {
    for (int m=0; m < unknowns; m++) {
      const int i_Q = i*unknowns  + m;
      validatePicardLoopHullAndCorrectorResult(Q,i_Q,m);
      buffer << std::fixed << std::showpoint << std::setprecision(6) << Q[i_Q] << " ";
    }
    buffer << "\n";
  }
  if ( verbose ) {
    std::cout << buffer.str() << std::flush; 
  }
  
  ::exahype2::aderdg::corrector_addCellContributions_loop_AoS(
    [&] (
      const double * const __restrict__           Q,
      const tarch::la::Vector<Dimensions,double>& x,
      double                                      t,
      int                                         direction,
      double * __restrict__                       F
    ) -> void {
      flux(Q,x,t,direction,F);
    },
    [&] (
      const double * const __restrict__           Q,
      const tarch::la::Vector<Dimensions,double>& x,
      double                                      t,
      double * __restrict__                       S
    ) -> void {
    },
    [&] (
      const double * const __restrict__           Q,
      double * __restrict__                       dQ_or_deltaQ,
      const tarch::la::Vector<Dimensions,double>& x,
      double                                      t,
      int                                         direction,
      double * __restrict__                       BgradQ
    ) -> void {
    },
    U,                                      // UOut,                               
    Q,                                      // QIn, 
    QuadraturePoints,                       // nodes,
    QuadratureWeights,                      // weights,
    StiffnessOperator,                      // Kxi,
    DerivativeOperator,                     // dudx, 
    x,                                      // cellCentre,
    dx,                                     // dx,
    t,                                      // t,
    dt,                                     // dt,
    order,                                  // order,
    unknowns,                               // unknowns,
    auxiliaryVariables,                     // auxiliaryVariables,
    true /*callFlux*/,                      // callFlux,
    false /*callSource*/,                   // callSource,
    false /*callNonconservativeProduct*/
  );  // callNonconservativeProduct) 


  // TODO Using simple extrapolation leads to wrong result
  //::exahype2::aderdg::spaceTimePredictor_extrapolateInTime_loop_AoS(
  //  U,
  //  Q,
  //  BasisFunctionValuesRight,           // FRCoeff,
  //  order,
  //  unknowns,
  //  auxiliaryVariables);

  // print result
  buffer << "\nresult (U, post volume integral):\n\n";
  for (int i = 0; i < nodesPerCell; i++) {
    for (int m=0; m < unknowns; m++) {
      const int i_U = i*unknowns  + m;
      buffer << std::fixed << std::showpoint << std::setprecision(6) << U[i_U] << " ";
    }
    buffer << "\n";
  }
  if ( verbose ) {
    std::cout << buffer.str() << std::flush; 
  }
 
  // prepare QHull
  for (int face = 0; face < 2*Dimensions; face++) {
    for (int i = 0; i < 2*nodesPerCell; i++) {
      for (int m=0; m < unknowns; m++) {
        const int i_QHull = i*unknowns  + m;
        QHull[face][i_QHull] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

  ::exahype2::aderdg::spaceTimePredictor_extrapolate_loop_AoS(
    QHull,
    Q,
    BasisFunctionValuesLeft,            // FLCoeff,
    BasisFunctionValuesRight,           // FRCoeff,
    order,
    unknowns,
    auxiliaryVariables);
  
  // print result
  buffer << "\nresult (QHull):\n\n";
  for (int face = 0; face < 2*Dimensions; face++) {
    for (int i = 0; i < 2*nodesPerCell; i++) {
      for (int m=0; m < unknowns; m++) {
        const int i_QHull = i*unknowns  + m;
        const double value = QHull[face][i_QHull];
        buffer << std::fixed << std::showpoint << std::setprecision(6) << value << " ";
      }
      buffer << "\n";
    }
    buffer << "\n";
  }
  if ( verbose ) {
    std::cout << buffer.str() << std::flush; 
  }
  
  // extrapolate again with swapped face pointers (left<->right) to fill the neigbour blocks
  ::exahype2::aderdg::spaceTimePredictor_extrapolate_loop_AoS(
    QHull_faceSwap,
    Q,
    BasisFunctionValuesLeft,            // FLCoeff,
    BasisFunctionValuesRight,           // FRCoeff,
    order,
    unknowns,
    auxiliaryVariables);
  
  // print result 
  buffer << "\nresult (QHull, post filling neighbour arrays):\n\n";
  for (int face = 0; face < 2*Dimensions; face++) {
    for (int i = 0; i < 2*nodesPerCell; i++) {
      for (int m=0; m < unknowns; m++) {
        const int i_QHull = i*unknowns  + m;
        const double value = QHull[face][i_QHull];
        validatePicardLoopHullAndCorrectorResult(QHull[face],i_QHull,m);
        buffer << std::fixed << std::showpoint << std::setprecision(6) << value << " ";
      }
      buffer << "\n";
    }
    buffer << "\n";
  }
  if ( verbose ) {
    std::cout << buffer.str() << std::flush; 
  }

  // max eigenvalue
  ::exahype2::aderdg::riemann_maxAbsoluteEigenvalue_loop_AoS(
    maxEigenvalue,
    maxEigenvaluePerFace,
    QHull,
    QuadraturePoints,
    x,
    dx,
    t,
    dt,
    order,
    unknowns,
    auxiliaryVariables);
  
  buffer << "\nresult (maxEigenvaluePerFace):\n\n";
  for (int face = 0; face < 2*Dimensions; face++) {
    buffer << std::fixed << std::showpoint << std::setprecision(6) << maxEigenvaluePerFace[face] << " ";
  }
  if ( verbose ) {
    std::cout << buffer.str() << std::flush; 
  }
  // Rusanov test
  ::exahype2::aderdg::rusanovNonlinear_loop_AoS(
    flux,
    flux,
    [&](
      double * __restrict__                       Q,
      double * __restrict__                       dQ_or_deltaQ,
      const tarch::la::Vector<Dimensions,double>& x,
      double                                      t,
      int                                         direction,
      double * __restrict__                       BgradQ
    ) -> void {},
    [&](
      double * __restrict__                       Q,
      double * __restrict__                       dQ_or_deltaQ,
      const tarch::la::Vector<Dimensions,double>& x,
      double                                      t,
      int                                         direction,
      double * __restrict__                       BgradQ
    ) -> void {},
    riemannResult,
    QHull, 
    maxEigenvaluePerFace,
    QuadraturePoints, 
    QuadratureWeights, 
    x,
    dx,
    t,
    dt,
    order,
    unknowns,
    auxiliaryVariables,
    atBoundary,
    true /*callFlux*/,
    false/*callNonconservativeProduct*/);

  // print result
  buffer << "\nresult (riemannResult):\n\n";
  for (int face = 0; face < Dimensions*2; face++) {
    for (int i = 0; i < nodesPerFace; i++) {
      for (int m=0; m < unknowns; m++) {
        const int i_QFace = (face*nodesPerFace+i)*unknowns + m;
        buffer << std::fixed << std::showpoint << std::setprecision(6) << riemannResult[i_QFace] << " ";
      }
      buffer << "\n";
    }
    buffer << "\n";
  }
  if ( verbose ) {
    std::cout << buffer.str() << std::flush; 
  }

  double maxEigenvalueInCell = ::exahype2::aderdg::corrector_addRiemannContributions_loop_AoS(
    [&] (
      double * __restrict__                       Q,
      const tarch::la::Vector<Dimensions,double>& x,
      double                                      t
    ) -> void {
    },
    maxEigenvalue,
    U,
    riemannResult,
    QuadraturePoints,
    QuadratureWeights,
    BasisFunctionValuesLeft, // only left
    x,
    dx,
    t,
    dt,
    order,
    unknowns,
    auxiliaryVariables,
    false); // callMaxEigenvalue
  
  buffer << "\nresult (U, post Riemann):\n\n";
  for (int i = 0; i < nodesPerCell; i++) {
    for (int m=0; m < unknowns; m++) {
      const int i_U = i*unknowns  + m;
      validatePicardLoopHullAndCorrectorResult(U,i_U,m);
      buffer << std::fixed << std::showpoint << std::setprecision(6) << U[i_U] << " ";
    }
    buffer << "\n";
  }
  if ( verbose ) {
    std::cout << buffer.str() << std::flush; 
  }
}

void exahype2::aderdg::tests::ADERDGTest::run() {
  testMethod (testAdvection)
  testMethod (testEuler)
  testMethod (testInterpolate)
}

void exahype2::aderdg::tests::ADERDGTest::testInterpolate() {
  std::ostringstream buffer;
  const bool verbose = false;
  
  constexpr int unknowns           = 5;
  constexpr int auxiliaryVariables = 0;
  constexpr int order              = 3; // order must match nodes, weights etc.

  constexpr int strideQ       = unknowns+auxiliaryVariables;
  constexpr int nodesPerAxis  = (order+1);
  int nodesPerFace  = nodesPerAxis; 
  if ( Dimensions == 3 ) {
    nodesPerFace *= nodesPerAxis;
  }
  const int nodesPerCell = nodesPerFace*nodesPerAxis;

  if ( verbose ) {
    std::cout << "# INTERPOLATION:\n" << std::endl;
  }
  
  // In-/Outputs:
  buffer << "\ninput (U):\n\n";
  double U[nodesPerCell*strideQ];
  for (int i = 0; i < nodesPerCell; i++) {
    for (int m=0; m < unknowns; m++) {
      const int i_U = i*strideQ + m;
      U[i_U] = testEuler_UIn[m];
      buffer << std::fixed << std::showpoint << std::setprecision(6) << U[i_U] << " ";
    }
    for (int m=unknowns; m < strideQ; m++) {
      const int i_U = i*strideQ + m;
      U[i_U] = 0.0;
      buffer << std::fixed << std::showpoint << std::setprecision(6) << U[i_U] << " ";
    }
    buffer << "\n";
  }
  if ( verbose ) {
    std::cout << buffer.str() << std::flush; 
  }
  
  const tarch::la::Vector<Dimensions,double> refX(0.5); // must be in [0,1]^d

  double pointwiseQOut[strideQ]; // unknowns + auxiliaryVariables
  ::exahype2::aderdg::interpolate_AoS(
    U,
    QuadraturePoints,
    BarycentricWeights,
    refX,
    nodesPerAxis,
    strideQ,
    pointwiseQOut);

  buffer << "\noutput (pointwiseQOut):\n\n";
  for (int m=0; m < unknowns; m++) {
    const double value = pointwiseQOut[m];
      const double eps = 1.0e-6;
      validateNumericalEqualsWithEpsWithParams1(
        value, testEuler_UIn[m], eps, m); // constant solution
    buffer << std::fixed << std::showpoint << std::setprecision(6) << value << " ";
  }
  buffer << "\n";
  if ( verbose ) {
    std::cout << buffer.str() << std::flush; 
  }
}

void exahype2::aderdg::tests::ADERDGTest::testAdvection() {
  // geometry
  const tarch::la::Vector<Dimensions, double> x({0.0, 0.0});
  const double dx = 1;
  const double t  = 0.0;
  const double dt = 0.001;

  const int unknowns           = 4;
  const int auxiliaryVariables = 0;
  const int order              = 3; // order must match nodes, weights etc.

  const bool verbose = false;

  if ( verbose ) {
    std::cout << "# LINEAR ADVECTION:\n" << std::endl;
  }
  runADERDGStep(
    [&] (      
      const double * const __restrict__           Q,
      const tarch::la::Vector<Dimensions,double>& x,
      double                                      t,
      int                                         direction,
      double * __restrict__                       F
    ) -> void {
      const double irho = 1./Q[0];
      const double velocity = irho*Q[direction+1];
      F[0] = velocity*Q[0];
      F[1] = 0.0;
      F[2] = 0.0;
      F[3] = 0.0;
    },
    [&] (
      const double * const __restrict__           Q,
      const tarch::la::Vector<Dimensions,double>& x,
      double                                      t,
      const int                                   direction
    ) -> double {
      const double irho = 1./Q[0];
      const double velocity = irho*Q[direction+1];
      return std::abs(velocity);
    },
    [&] (
      const double* Q,
      const int     i_Q,
      const int     m 
    ) -> void { 
      const double eps = 1.0e-6;
      validateNumericalEqualsWithEpsWithParams1(
        Q[i_Q], testAdvection_UIn[m], eps, i_Q); // constant solution
    },
    x,dx,t,dt,unknowns,auxiliaryVariables,order,
    testAdvection_UIn, 
    verbose
  );
}

void exahype2::aderdg::tests::ADERDGTest::testEuler() {
  // geometry
  const tarch::la::Vector<Dimensions, double> x({0.0, 0.0});
  const double dx = 5e-02;
  const double t  = 0.0;
  const double dt = 1.686854344081342E-003;

  const int unknowns           = 5;
  const int auxiliaryVariables = 0;
  const int order              = 3; // order must match nodes, weights etc.

  const bool verbose = false;

  if ( verbose ) { 
    std::cout << "\n# COMPRESSIBLE EULER:\n" << std::endl;
  }
  runADERDGStep(
    [&] (      
      const double * const __restrict__           Q,
      const tarch::la::Vector<Dimensions,double>& x,
      double                                      t,
      int                                         direction,
      double * __restrict__                       F
    ) -> void {
      constexpr double gamma = 1.4;
      const double irho = 1./Q[0];
      const double p = (gamma-1) * (Q[4] - 0.5*irho*(Q[1]*Q[1]+Q[2]*Q[2]));
    
      const double velocity = irho*Q[direction+1];
      F[0] = velocity*Q[0];
      F[1] = velocity*Q[1];
      F[2] = velocity*Q[2];
      F[3] = velocity*Q[3];
      F[direction+1] += p;
      F[4] = velocity*(Q[4]+p);  
    },
    [&] (
      const double * const __restrict__           Q,
      const tarch::la::Vector<Dimensions,double>& x,
      double                                      t,
      const int                                   direction
    ) -> double {
      constexpr double gamma = 1.4;
      const double irho = 1./Q[0];
      #if Dimensions==3
      const double p = (gamma-1) * (Q[4] - 0.5*irho*(Q[1]*Q[1]+Q[2]*Q[2]+Q[3]*Q[3]));
      #else
      const double p = (gamma-1) * (Q[4] - 0.5*irho*(Q[1]*Q[1]+Q[2]*Q[2]));
      #endif
    
      const double u_n = Q[direction + 1] * irho;
      const double c   = std::sqrt(gamma * p * irho);
    
      return std::max( std::abs(u_n - c), std::abs(u_n + c) );
    },
    [&] (
      const double* Q,
      const int     i_Q,
      const int     m
    ) -> void { 
      const double eps = 1.0e-5;
      validateNumericalEqualsWithEpsWithParams1(
        Q[i_Q], testEuler_QOut[m], eps, i_Q);
    },
    x,dx,t,dt,unknowns,auxiliaryVariables,order,
    testEuler_UIn, 
    verbose
  );
}
