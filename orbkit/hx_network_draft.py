'''Drafts for ZIBAmira network files'''

#: Draft for a standard ZIBAmira network
hx_network= '''
# Amira Script
remove -all
remove {FILENAME.am} {Result} {FILENAME.cmap} {Arithmetic} {Voltex}

# Create viewers
viewer setVertical 0

viewer 0 setBackgroundMode 0
viewer 0 setBackgroundColor 1 1 1
viewer 0 setBackgroundColor2 1 1 1
viewer 0 setTransparencyType 5
viewer 0 setAutoRedraw 0
viewer 0 show
mainWindow show

set hideNewModules 0
[ load ${SCRIPTDIR}/FILENAME.am ] setLabel {FILENAME.am}
{FILENAME.am} setIconPosition 20 10
FILENAME.am fire
FILENAME.am setViewerMask 16383

set hideNewModules 0
[ load ${SCRIPTDIR}/FILENAME.cmap ] setLabel {FILENAME.cmap}
{FILENAME.cmap} setIconPosition 20 71
{FILENAME.cmap} fire
FILENAME.cmap fire
FILENAME.cmap setViewerMask 16383

set hideNewModules 0
create {HxArithmetic} {Arithmetic}
{Arithmetic} setIconPosition 160 10
{Arithmetic} {inputA} connect FILENAME.am
{Arithmetic} fire
{Arithmetic} {resultChannels} setIndex 0 0
{Arithmetic} {nValues} setMinMax 0 -2147483648 2147483648
{Arithmetic} {nValues} setValue 0 1
{Arithmetic} {expr0} setState A
{Arithmetic} {resultType} setValue 1
{Arithmetic} {resultLocation} setValue 0
{Arithmetic} {resolution} setMinMax 0 1 100000
{Arithmetic} {resolution} setValue 0 100
{Arithmetic} {resolution} setMinMax 1 1 100000
{Arithmetic} {resolution} setValue 1 100
{Arithmetic} {resolution} setMinMax 2 1 100000
{Arithmetic} {resolution} setValue 2 100
{Arithmetic} {minBox} setMinMax 0 -3.40282346638529e+38 3.40282346638529e+38
{Arithmetic} {minBox} setValue 0 -1
{Arithmetic} {minBox} setMinMax 1 -3.40282346638529e+38 3.40282346638529e+38
{Arithmetic} {minBox} setValue 1 -1
{Arithmetic} {minBox} setMinMax 2 -3.40282346638529e+38 3.40282346638529e+38
{Arithmetic} {minBox} setValue 2 -1
{Arithmetic} {maxBox} setMinMax 0 -3.40282346638529e+38 3.40282346638529e+38
{Arithmetic} {maxBox} setValue 0 1
{Arithmetic} {maxBox} setMinMax 1 -3.40282346638529e+38 3.40282346638529e+38
{Arithmetic} {maxBox} setValue 1 1
{Arithmetic} {maxBox} setMinMax 2 -3.40282346638529e+38 3.40282346638529e+38
{Arithmetic} {maxBox} setValue 2 1
{Arithmetic} {options} setValue 0 0
{Arithmetic} {doIt} snap 0 1
Arithmetic fire
Arithmetic setViewerMask 16383
Arithmetic setPickable 1

set hideNewModules 0
[ {Arithmetic} create result ] setLabel {Result}
{Result} setIconPosition 21 41
{Result} {master} connect {Arithmetic} result
Result fire
Result setViewerMask 16383

set hideNewModules 0
create {HxVoltex} {Voltex}
{Voltex} setIconPosition 470 41
{Voltex} {data} connect Result
{Voltex} {colormap} setDefaultColor 1 0.8 0.5
{Voltex} {colormap} setDefaultAlpha 0.500000
{Voltex} {colormap} setLocalRange 0
{Voltex} {colormap} connect FILENAME.cmap
{Voltex} fire
{Voltex} {options} setValue 0 0
{Voltex} {options} setValue 1 1
{Voltex} {range} setMinMax 0 -3.40282346638529e+38 3.40282346638529e+38
{Voltex} {range} setValue 0 0
{Voltex} {range} setMinMax 1 -3.40282346638529e+38 3.40282346638529e+38
{Voltex} {range} setValue 1 255
{Voltex} {lookup} setValue 2
{Voltex} {gamma} setMinMax 1 8
{Voltex} {gamma} setButtons 0
{Voltex} {gamma} setIncrement 0.466667
{Voltex} {gamma} setValue 3
{Voltex} {gamma} setSubMinMax 1 8
{Voltex} {alphaScale} setMinMax 0 1
{Voltex} {alphaScale} setButtons 0
{Voltex} {alphaScale} setIncrement 0.1
{Voltex} {alphaScale} setValue 1
{Voltex} {alphaScale} setSubMinMax 0 1
{Voltex} {texture2D3D} setValue 1
{Voltex} {slices} setMinMax 0 5120
{Voltex} {slices} setButtons 1
{Voltex} {slices} setIncrement 1
{Voltex} {slices} setValue 512
{Voltex} {slices} setSubMinMax 0 512
{Voltex} {downsample} setMinMax 0 1 100
{Voltex} {downsample} setValue 0 1
{Voltex} {downsample} setMinMax 1 1 100
{Voltex} {downsample} setValue 1 1
{Voltex} {downsample} setMinMax 2 1 100
{Voltex} {downsample} setValue 2 1
{Voltex} {doIt} snap 0 1
{Voltex} doIt hit
Voltex fire
Voltex setViewerMask 16383
Voltex select
Voltex setPickable 1

set hideNewModules 0


viewer 0 setCameraOrientation -0.869144 0.48907 -0.0734764 2.86511
viewer 0 setCameraPosition 15.0752 9.71314 -55.4207
viewer 0 setCameraFocalDistance 58.25
viewer 0 setCameraNearDistance 39.1308
viewer 0 setCameraFarDistance 77.4073
viewer 0 setCameraType orthographic
viewer 0 setCameraHeight 48.1395
viewer 0 setAutoRedraw 1
viewer 0 redraw
'''
