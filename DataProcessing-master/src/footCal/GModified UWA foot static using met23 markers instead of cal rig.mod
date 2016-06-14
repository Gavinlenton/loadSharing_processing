{********************************************************************************************}
{* Calculate ABD/ADD angle of foot for UWA foot CS definition in the static trial where,
   in place of UWA foot rig, LMT23/RMT23 markers are used and the heel markers (LCAL/RCAL)
   are aligned along either of the two horizontal axes of the global coordinate system
   (eg either edge of L-frame force plate) *}
{* v2.5 Tim Wrigley, University of Melbourne 13/1/06 *}
{********************************************************************************************}

debug=1 {* set to 0 for runtime *}

{*
Marker configuration with respect to global CS, looking down on floor plane:


		                     ^ global axis (X, Y or Z)
		                     |
		   LMT23   RMT23     |
		    *        *       |
		                     |
		                     |
		                     |
		                     |
		                     |
		    LCAL  RCAL       |
		<-----*----*---------+ global origin
 global axis (X, Y or Z)


Will also work if LCAL/RCAL aligned with other horizontal axis.
Does not assume any particular coordinate system (eg Vicon or ISB).
*}


{*********** Start of macro section **************}

Macro DrawGlobal(ScaleFactor)
{* draws global coordinate system *}
{* Oglobal, Xglobal, Yglobal, Zglobal must be in MKR file for display *}
Oglobal={0,0,0}
Xglobal=ScaleFactor*{1,0,0}
Yglobal=ScaleFactor*{0,1,0}
Zglobal=ScaleFactor*{0,0,1}
OUTPUT(Xglobal,Yglobal,Zglobal,Oglobal)
EndMacro

Macro DrawCS(segment,scale)
{* draws segment coordinate system *}
{* O#segment, X#segment, Y#segment, Z#segment must be in MKR file for display *}
O#segment={0,0,0}*segment
X#segment=O#segment+scale*1(segment) 
Y#segment=O#segment+scale*2(segment) 
Z#segment=O#segment+scale*3(segment) 
Output(O#segment,X#segment,Y#segment,Z#segment)
EndMacro

macro GetMaxVecCpt(vec,cpt,value) {* get the largest component of vector *}
value=vec(1)
cpt=1
if vec(2) > value
	value=vec(2)
	cpt=2
endif
if vec(3) > value
	value=vec(3)
	cpt=3
endif
endmacro

macro GetMinVecCpt(vec,cpt,value) {* get the smallest component of vector *}
value=vec(1)
cpt=1
if vec(2) < value
	value=vec(2)
	cpt=2
endif
if vec(3) < value
	value=vec(3)
	cpt=3
endif
endmacro

macro SetVecCpt(vec,cpt,value) {* set cpt of vec to value *}
{* macro required since BB does not allow "vec(cpt) = value" *}
if cpt == 1
	vec={value,vec(2),vec(3)}
elsif cpt == 2
	vec={vec(1),value,vec(3)}
else
	vec={vec(1),vec(2),value}
endif
endmacro

macro GetVecCpt(vec,cpt,value) {* get value of cpt of vec *}
{* macro required since BB does not allow "variable = vec(cpt)" where cpt is a variable *}
if cpt == 1
	value = vec(1)
elsif cpt == 2
	value = vec(2)
else
	value = vec(3)
endif
endmacro

macro GetCardanAngles(angles, childsegment, parentsegment, order, bodyfixed)
{* bodyfixed=0 => space-fixed sequence (BB default); bodyfixed=1 => body-fixed sequence (Kane terminology) *}
{* Gets Cardan angles for either space- or body-fixed sequence, without resorting to inelegant 'trick' described in BB docs for body-fixed. *}
{* Also, BB only allows rotation order to be specified as numeric value (-3 to 3) or text token (xyz etc). It will not accept value or text token in 
   a variable; this macro allows order to be passed as a scalar variable. *}

if bodyfixed == 1 {* reverse the rotation sequence order for body-fixed *}
	if order == 1     {* xyz *}
		order=-3  {* ... zyx *}
	elsif order == 2  {* yzx *}
		order=-1  {* ... xzy *}
	elsif order == 3  {* zxy *}
		order=-2  {* ... yxz *}
	elsif order == -1 {* xzy *}
		order=2   {* ... yzx *}
	elsif order == -2 {* yxz *}
		order=3   {* ... zxy *}
	elsif order == -3 {* zyx *}
		order=1   {* ... xyz *}
	endif
endif
{* Get the angles *}
if order == 1
	angles = <childsegment,parentsegment,1>
elsif order == 2
	angles = <childsegment,parentsegment,2>
elsif order == 3
	angles = <childsegment,parentsegment,3>
elsif order == -1
	angles = <childsegment,parentsegment,-1>
elsif order == -2
	angles = <childsegment,parentsegment,-2>
elsif order == -3
	angles = <childsegment,parentsegment,-3>
endif

if bodyfixed == 1 {* reverse output order so angles are in originally specified body-fixed order *}
	angles = <angles(3), angles(2), angles(1)>
endif

endmacro


{**************** End of macro section ******************}


{**************** Algorithm section *********************}



{* Only runs if met23 markers present *}
if ExistAtAll(LMT23,RMT23)


{********** STEP 1: Find which global axis the heel markers are aligned along *}

{* Smallest vec cpt of mid-point between met markers will be vertical axis cpt.
   Assumes: foot is on the floor, ie near 0 plane for X, Y, or Z *}
{* MidMets=average((LMT23-RMT23)/2) *} {* average not allowed within IF block *}
$MidMets=RMT23+(LMT23-RMT23)/2
PARAM($MidMets) {* PARAM has implicit average, so inelegant way to get an AVERAGE done in an IF block *}
GetMinVecCpt( ABS($MidMets), VertAxisCpt, val) {* $ on MidMets so that it retrieves it from MP file *}
if debug PARAM( VertAxisCpt ) endif

{* find which are the 2 horizontal axes (axisA, axisB) *}
if VertAxisCpt == 1 {* X *}
	axisA = 2   {* Y *}
	axisB = 3   {* Z *}
elsif VertAxisCpt == 2
	axisA = 1
	axisB = 3
elsif VertAxisCpt == 3
	axisA = 1
	axisB = 2
endif

{* Create virtual heels at same height as L/RMT23 markers - ensures that foot CS will have vertical axis parallal to global vertical *}
VLCAL = LCAL
VRCAL = RCAL
GetVecCpt(LMT23,VertAxisCpt ,tmp1) {* height of LMT23 *}
GetVecCpt(RMT23,VertAxisCpt ,tmp2)
SetVecCpt(VLCAL ,VertAxisCpt , tmp1) {* set height of VLCAL to that of LMT23 *}
SetVecCpt(VRCAL ,VertAxisCpt , tmp2)
if debug param(VLCAL, VRCAL) endif

{* Longest distance of mid heel markers to horizontal axes is the axis along which the heels were aligned *}
{* Assumes: heel markers are aligned along (ie very close to) either of the two global horizontal axes *}
{* MidHeels = average((LCAL-RCAL)/2) *} {* average not allowed within IF block *}
MidHeels = (LCAL-RCAL)/2
MidHeelsX = ABS(MidHeels(1))
MidHeelsY = ABS(MidHeels(2))
MidHeelsZ = ABS(MidHeels(3))
$MidHeels = {MidHeelsX,MidHeelsY,MidHeelsZ}
PARAM($MidHeels) {* get the average done *}

{* BB doesn't allow this "MidHeels( VertAxisCpt ) = 1000000", so use this ... *}
SetVecCpt($MidHeels,VertAxisCpt, 1000000) {* set vert axis cpt to big num so min doesn't find it *}
GetMinVecCpt($MidHeels, SmallestAxisCpt, tmp )

Param(SmallestAxisCpt,axisA,axisB)
if SmallestAxisCpt == axisA
	BaselineAxis = axisB {* axisA cpt = distance FROM axisB (eg if X cpt smallest, then heels nearest to Y axis) *}
else
	BaselineAxis = axisA
endif
PARAM(BaselineAxis)
{********** STEP 2: Calculate the vector along the long axis of the foot, and the 'baseline' axis along which heels were aligned *}

LeftFootLine = LMT23-VLCAL {* use virtual heels (same height as LMT23) so line is horizontal *}
RightFootLine = RMT23-VRCAL
if BaselineAxis == 1 {* X *}
	Baseline = {1, 0, 0}
elsif BaselineAxis == 2 {* Y *}
	Baseline = {0, 1, 0}
else {* Z *}
	Baseline = {0, 0, 1}
endif
if debug param(LeftFootLine , RightFootLine) endif
PARAM(BaselineAxis)

{********** STEP 3: Calculate the foot ABD/ADD angle, between the long foot axis and the axis along which the heels are aligned *}

{* Foot segment definition:
   Origin at L/RCAL
   1st (x) axis through foot long axis, but horizontal
   2nd (z) axis is cross of global axis nearest heels onto 1st axis
   3rd (y) axis is cross of 2nd onto 1st *}
{* If RCAL closer to global origin on baseline axis, need to reverse direction of baseline vector to give vert axis-up for this cross
   (if LCAL closer to origin on baseline axis, then cross already gives vert axis-up) *}
GetVecCpt(RCAL,BaselineAxis,tmp1)
GetVecCpt(LCAL,BaselineAxis,tmp2)
if debug param(tmp1,tmp2) endif
if tmp1 < tmp2
	Baseline = -Baseline
endif
if debug param(Baseline) endif

LFOOTANGCS = [LCAL, LeftFootLine, Baseline, yzx]
RFOOTANGCS = [RCAL, RightFootLine, Baseline, yzx]

{* To get foot ABD/ADD projected angle, use the equivalence between projected angle and one angle
   in a particular Cardan space-fixed sequence (Crawford et al 96 Hum Movt Sci):
Pxj = Ryzx(3) ie projected angle about x (onto yz plane) = x angle, last in Ryzx Cardan space-fixed sequence
Pxk = Rzyx(3)
Pyi = Rxzy(3)
Pyk = Rzxy(3)
Pzi = Rxyz(3)
Pzj = Ryxz(3)
(i j k refers to the orthogonal segment axes (2 in each case) to which each projected angle can be referenced)
*}
 
{* Projected angle also could be done with dot product or cross product of two vectors (with zeroing of vector cpts in the projection plane,
   but easier to get BB to do the work via Cardan angles ... (?) *}

{* Set token for Cardan space-fixed sequence that will give third angle equivalent to rotation angle about vertical axis (projected onto horiz plane):
	1. Rotation about horizontal axis (neutral long axis of foot) opposite to baseline axis
	2. Rotation about baseline axis
	3. Rotation about vertical axis (this is equivalent to projected angle onto horizontal plane) *}
if VertAxisCpt == 1 AND BaselineAxis == 2
	order = -3 {* zyx = 321 *}
elsif VertAxisCpt == 1 AND BaselineAxis == 3
	order = 2 {* yzx = 231 *}
elsif VertAxisCpt == 2 AND BaselineAxis == 1
	order = 3 {* zxy = 312 *}
elsif VertAxisCpt == 2 AND BaselineAxis == 3
	order = -1 {* xzy = 132 *}
elsif VertAxisCpt == 3 AND BaselineAxis == 1
	order = -2 {* yxz = 213 *}
elsif VertAxisCpt == 3 AND BaselineAxis == 2
	order = 1 {* xyz = 123 *}
endif
if debug param(order) endif


GetCardanAngles(LFOOTROT,LFOOTANGCS,1,order,0) {* parent = 1 => global CS *}
GetCardanAngles(RFOOTROT,RFOOTANGCS,1,order,0)

{* Projected angle = last angle in space-fixed Cardan sequence (rotation about vertical axis) *}
{* Positive according to RH screw about vertical axis *}
$LFootAbduction=LFOOTROT(3) {* ABD positive *}
$RFootAbduction=-RFOOTROT(3)

{************* Output section ************}

PARAM($LFootAbduction) {* '$' to make these appear as Subject Measurements in WS *}
PARAM($RFootAbduction)


{************* Drawing section ************}

DrawGlobal(1000) {* Vicon Z-up global CS *}
DrawCS( LFOOTANGCS, 500)
DrawCS( RFOOTANGCS, 500)


endif {* end If ExistAtAll(LMT23,RMT23) *}
