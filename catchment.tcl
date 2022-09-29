#!/usr/bin/env wish
# ===========================================================================
# 
# catchment.tcl --
# plot catchment area
# 
# This file is part of the following software:
# 
#    - the low-level C++ template SIMD library
#    - the SIMD implementation of the MinWarping and the 2D-Warping methods 
#      for local visual homing.
# 
# The software is provided based on the accompanying license agreement
# in the file LICENSE or LICENSE.doc. The software is provided "as is"
# without any warranty by the licensor and without any liability of the
# licensor, and the software may not be distributed by the licensee; see
# the license agreement for details.
# 
# (C) Ralf Möller
#     Computer Engineering
#     Faculty of Technology
#     Bielefeld University
#     www.ti.uni-bielefeld.de
# 
# ===========================================================================

# code adapted from:
# ~/PROG/applications/Homing-6-0/src:
#   plotCatchmentArea.tcl, plot.tcl, homing_init.tcl etc.
# ~/PROG/applications/Warping-6-0/src:
#   PlotCatchmentArea.tcl

# ===========================================================================
# assert
# ===========================================================================

# from https://wiki.tcl-lang.org/page/Assertions
proc assert condition {
    if {![uplevel 1 expr $condition]} {
        return -code error "assertion failed: $condition"
    }
}

# ===========================================================================
# flags
# ===========================================================================

set pcaShowGoalPoint 0x01
set pcaShowMovementLine 0x02
set pcaMarkSuccessfulPoints 0x04
set pcaPlotWithGreyBoxes 0x08
set pcaShowHomeVector 0x10
set pcaDemoMode 0x20
set pcaPrintStateArray 0x40

# ===========================================================================
# canvas parameter
# ===========================================================================

# not checked which are needed
set cvBox 40
set cvBox2 [expr $cvBox/2]
set cvArrowWidth 3
set cvArrowLen [expr $cvBox2]
set cvArrowShape {8 8 4}
set cvLineWidth 2
set cvCircleWidth 2
set cvCompassLen 20
set cvCompassWidth 5
# for drawVector2
set cvArrowWidth1 3
set cvArrowWidth2 [expr $cvArrowWidth1 + 3]
set cvArrowLen2 [expr $cvBox2]
set cvArrowLen1 [expr $cvArrowLen2 - 5]
set cvArrowShape1 {6 6 3}
set cvArrowShape2 {11 11 5}
set cvArrowBackOffset 2

# ===========================================================================
# catchment plot parameter
# ===========================================================================

set homeVectorColor blue3
set markerColor green2
set movementColor red2
set demoModeTimeDelay 200

set distThresh 0.6
# set goalMode "neighborhood"
set goalMode "point"
# set flags [expr $pcaPlotWithGreyBoxes | $pcaShowHomeVector]
set flags [expr $pcaShowGoalPoint | $pcaShowMovementLine | \
	       $pcaMarkSuccessfulPoints | $pcaShowHomeVector | \
	       $pcaDemoMode]

# ===========================================================================
# pointList functions
# ===========================================================================

proc cmp {e1 e2} {
    set c1 [lindex $e1 2]
    set c2 [lindex $e2 2]
    return [expr ($c1 < $c2) ? -1 : ($c1 == $c2) ? 0 : 1]
}

# make a point list in x-y order, list contains grid coordinates and
# grid distance from xss, yss
proc makePointList {xss yss} {
    global xmin ymin xmax ymax
    set pointList {}
    for {set x $xmin} {$x <= $xmax} {incr x} {
	for {set y $ymin} {$y <= $ymax} {incr y} {
	    # exclude xss,yss
	    if {($x != $xss) || ($y != $yss)} {
		set d [expr hypot(double($x - $xss), double($y - $yss))]
		lappend pointList "$x $y $d"
	    }
	}
    }
    return $pointList
}

# make a sorted list of images, starting with the ones distant to ss
proc makeSortedPointListReverse {xss yss} {
    set pointList [makePointList $xss $yss]
    set pointList [lsort -decreasing -command cmp $pointList]
    return $pointList
}

# ===========================================================================
# plot functions
# ===========================================================================

#-----------------------------------------------------------
# create canvas
#-----------------------------------------------------------

proc createCanvas {w title} {
    global cvBox xmax ymax xmin ymin

    set windowExists [winfo exists $w]
    if {$windowExists} {
	destroy $w
    }
    toplevel $w
    wm title $w $title
    canvas ${w}.cv -width [expr $cvBox * ($xmax+1-$xmin)] \
	-height [expr $cvBox * ($ymax+1-$ymin)] \
	-background white
    pack ${w}.cv -expand yes -fill both
    update
}

#-----------------------------------------------------------
# draw box in canvas (for marking home position)
#-----------------------------------------------------------

proc drawBox {w x y color {tags "box"}} {
    global cvBox cvLineWidth xmin ymin ymax

    # box coordinates
    set cvXtl [expr ($x - $xmin) * $cvBox]
    #  9. Dec 21 (rm)
    # set cvYtl [expr ($y - $ymin) * $cvBox]
    set cvYtl [expr ($ymax - $y) * $cvBox]
    set cvXbr [expr $cvXtl + $cvBox]
    set cvYbr [expr $cvYtl + $cvBox]
    # draw box
    ${w}.cv create rectangle $cvXtl $cvYtl $cvXbr $cvYbr \
	-outline $color -width $cvLineWidth -fill "" -tags $tags
    update
}

#-----------------------------------------------------------
# draw filled box in canvas
#-----------------------------------------------------------

proc fillBox {w x y color {tags "box"}} {
    global cvBox cvLineWidth xmin ymin ymax

    # box coordinates
    set cvXtl [expr ($x - $xmin) * $cvBox]
    #  9. Dec 21 (rm)
    # set cvYtl [expr ($y - $ymin) * $cvBox]
    set cvYtl [expr ($ymax - $y) * $cvBox]
    set cvXbr [expr $cvXtl + $cvBox]
    set cvYbr [expr $cvYtl + $cvBox]
    # draw box
    ${w}.cv create rectangle $cvXtl $cvYtl $cvXbr $cvYbr \
	-outline $color -width 1 -fill $color -tags $tags
    update
}

#-----------------------------------------------------------
# draw home vector in canvas
#-----------------------------------------------------------

# home vector is given in image coordinate system
#
# y
# ^
# |
# +--> x
proc drawVector {w x y xHome yHome color \
		     {arrowLen 0} {arrowWidth 0} {arrowShape {}} \
		     {tags "vec"}} {
    global cvBox cvBox2 cvArrowLen cvArrowShape \
	cvArrowWidth xmin ymin ymax

    if {$arrowLen == 0} {
	set arrowLen $cvArrowLen
    }
    if {$arrowWidth == 0} {
	set arrowWidth $cvArrowWidth
    }
    if {[llength $arrowShape] == 0} {
	set arrowShape $cvArrowShape
    }
    # home vector
    set cvX1 [expr ($x - $xmin) * $cvBox + $cvBox2]
    #  9. Dec 21 (rm)
    # set cvY1 [expr ($y - $ymin) * $cvBox + $cvBox2]
    set cvY1 [expr ($ymax - $y) * $cvBox + $cvBox2]
    #  5. Feb 06 (rm): that was -/- before
    set cvX2 [expr $cvX1 + $xHome * $arrowLen]
    set cvY2 [expr $cvY1 - $yHome * $arrowLen]
    ${w}.cv create line \
	$cvX1 $cvY1 \
	$cvX2 $cvY2 \
	-fill $color -width $arrowWidth \
	-arrow last -arrowshape $arrowShape -tags $tags
    update
}

#-----------------------------------------------------------
# draw line in canvas
#-----------------------------------------------------------

proc drawLine {w x1 y1 x2 y2 color {width 0} {tags "line"}} {
    global cvBox cvBox2 cvLineWidth xmin ymin ymax

    if {$width == 0} {
	set width $cvLineWidth
    }
    # transform to canvas coordinates
    set cvX1 [expr ($x1 - $xmin) * $cvBox + $cvBox2]
    #  9. Dec 21 (rm)
    # set cvY1 [expr ($y1 - $ymin) * $cvBox + $cvBox2]
    set cvY1 [expr ($ymax - $y1) * $cvBox + $cvBox2]
    set cvX2 [expr ($x2 - $xmin) * $cvBox + $cvBox2]
    #  9. Dec 21 (rm)
    # set cvY2 [expr ($y2 - $ymin) * $cvBox + $cvBox2]
    set cvY2 [expr ($ymax - $y2) * $cvBox + $cvBox2]
    ${w}.cv create line $cvX1 $cvY1 $cvX2 $cvY2 \
	-fill $color -width $width \
	-tags $tags
    update
}

#-----------------------------------------------------------
# draw circle in canvas
#-----------------------------------------------------------

proc drawCircle {w x y r outline fill {tags "circle"}} {
    global cvBox cvBox2 cvCircleWidth xmin ymin ymax

    # transform to canvas coordinates
    set cvR [expr $r * $cvBox]
    set cvX1 [expr ($x - $xmin) * $cvBox + $cvBox2 - $cvR]
    #  9. Dec 21 (rm)
    # set cvY1 [expr ($y - $ymin) * $cvBox + $cvBox2 - $cvR]
    set cvY1 [expr ($ymax - $y) * $cvBox + $cvBox2 - $cvR]
    set cvX2 [expr $cvX1 + 2 * $cvR]
    set cvY2 [expr $cvY1 + 2 * $cvR]
    ${w}.cv create oval $cvX1 $cvY1 $cvX2 $cvY2 \
	-fill $fill -outline $outline -width $cvCircleWidth \
	-tags $tags
    update
}

#-----------------------------------------------------------
# finish plot
#-----------------------------------------------------------

proc finishPlot {w} {
    # (why is this update necessary?)
    update
    ${w}.cv raise "box"
    ${w}.cv raise "vec"
    ${w}.cv lower "line"
    update
}

#-----------------------------------------------------------
# clearCanvas
#-----------------------------------------------------------

proc clearCanvas {w {tag all}} {
    ${w}.cv delete $tag
    update
}

#-----------------------------------------------------------
# savePlot
#-----------------------------------------------------------

proc savePlot {w filename} {
    puts "saving plot $w to $filename"
    ${w}.cv postscript -colormode color \
	-file $filename
}

# ===========================================================================
# closestNeighbor
# ===========================================================================

# determines closest location among the 8 neighbors in the grid and the
# distance (smallest perp. distance)
# x y is the integer grid coordinate of the current location
# hx hy is the home vector (unit length!)
# result: xn yn dist (new integer grid cordinate, minimal distance)
proc closestNeighbor {x y hx hy} {
    global xmin xmax ymin ymax
    # search range (with border treatment)
    set xl [expr $x - 1]
    set xl [expr $xl < $xmin ? $xmin : $xl]
    set xr [expr $x + 1]
    set xr [expr $xr > $xmax ? $xmax : $xr]
    set yl [expr $y - 1]
    set yl [expr $yl < $ymin ? $ymin : $yl]
    set yr [expr $y + 1]
    set yr [expr $yr > $ymax ? $ymax : $yr]
    # loop through all neighbors
    set minDist 1E10
    set xnMin $x
    set ynMin $y
    for {set xn $xl} {$xn <= $xr} {incr xn} {
	for {set yn $yl} {$yn <= $yr} {incr yn} {
	    # exclude center itself
	    if {($xn != $x) || ($yn != $y)} {
		# difference vector
		set dx [expr $xn - $x]
		set dy [expr $yn - $y]
		# squared len
		set d2 [expr $dx * $dx + $dy * $dy]
		# projection of difference vector on line through home vector
		set p [expr $hx * $dx + $hy * $dy]
		# only positive projections are considered
		if {$p >= 0.0} {
		    # squared len of projection
		    set p2 [expr $p * $p]
		    #puts "$x $y $hx $hy [expr hypot($hx,$hy)] $xn $yn $d2 $p2"
		    # perpendicular distance on line through home vector
		    set dist [expr sqrt($d2 - $p2)]
		    # search for mininmal distance
		    if {$dist < $minDist} {
			set minDist $dist
			set xnMin $xn
			set ynMin $yn
		    }
		}
	    }
	}
    }
    return [list $xnMin $ynMin $minDist]
}

# ===========================================================================
# getNeighborhood3x3
# ===========================================================================

# get coordinate list of neighborhood (database margins are considered)
proc getNeighborhood3x3 {xc yc} {
    global xmin xmax ymin ymax
    set nbh {}
    for {set i -1} {$i <= 1} {incr i} {
	set x [expr $xc + $i]
	if {($x >= $xmin) && ($x <= $xmax)} {
	    for {set j -1} {$j <= 1} {incr j} {
		set y [expr $yc + $j]
		if {($y >= $ymin) && ($y <= $ymax)} {
		    lappend nbh [list $x $y]
		}
	    }
	}
    }
    return $nbh
}

# ===========================================================================
# plotCatchmentArea
# ===========================================================================

# plot catchment area in canvas
# goalMode: "point" or "neighborhood"
proc plotCatchmentArea {canvas xss yss distThresh goalMode flags} {
    global xmax xmin ymax ymin homeVectors
    global pcaShowGoalPoint pcaShowMovementLine pcaMarkSuccessfulPoints \
	pcaPlotWithGreyBoxes pcaShowHomeVector pcaDemoMode \
	pcaPrintStateArray
    global homeVectorColor markerColor movementColor demoModeTimeDelay

    set stateChar [list U T F S G]
    # number of steps in each run: 2 * diagonal length
    set maxSteps [expr 2.0 * hypot($xmax - $xmin, $ymax - $ymin)] 
    # visualize goals
    if {$flags & $pcaShowGoalPoint} {
	# draw snapshot position
	drawBox $canvas $xss $yss black
	drawCircle $canvas $xss $yss 0.05 black black
    }
    # different states of nodes
    set undecided 0
    set timeout 1
    set failure 2
    set success 3
    set goal 4
    # set state of each node to undecided
    array unset state
    for {set x $xmin} {$x <= $xmax} {incr x} {
	for {set y $ymin} {$y <= $ymax} {incr y} {
	    set state($x,$y) "$undecided"
	}
    }
    if {$goalMode == "neighborhood"} {
	# mark 3x3 neighborhood as goal
	foreach nb [getNeighborhood3x3 $xss $yss] {
	    lassign $nb xnb ynb
	    set state($xnb,$ynb) "$goal"
	}
    } elseif {$goalMode == "point"} {
	# mark only goal point as goal
	set state($xss,$yss) "$goal"
    } else {
	puts stderr "invalid goal mode $goalMode"
	exit
    }
    # draw all home vectors (but exclude goal position)
    for {set x $xmin} {$x <= $xmax} {incr x} {
	for {set y $ymin} {$y <= $ymax} {incr y} {
	    #  1. Jul 10 (rm): we exclude the goal positions
	    if {($x != $xss) || ($y != $yss)} {
		# draw home vector
		if {$flags & $pcaShowHomeVector} {
		    lassign $homeVectors($x,$y) xHomeW yHomeW
		    drawVector $canvas $x $y $xHomeW $yHomeW $homeVectorColor
		}
	    }
	}
    }
    # homing run
    # points with largest distance from snapshot location first
    set pointList [makeSortedPointListReverse $xss $yss]
    # puts $pointList
    # go through all starting points
    foreach point $pointList {
	# starting point
	lassign $point x y d
	set visitedPoints {}
	# if for-loop runs through, we also failed
	set marker "$timeout"
	# go through a predefined number of steps
	for {set step 0} {$step < $maxSteps} {incr step} {
	    # puts "$x $y"
	    # did we reach the goal yet?
	    if {$state($x,$y) == "$goal"} {
		# we have reached the goal
		if {$goalMode == "neighborhood"} {
		    # neighborhood point is added to list and marked below
		    lappend visitedPoints [list $x $y]
		}
		# mark all visited points as success later
		set marker "$success"
		# homing run is finished
		break
	    }
	    # we haven't reached the goal yet
	    if {$state($x,$y) == "$undecided"} {
		# append point to list
		lappend visitedPoints [list $x $y]
		# get precomputed home vector
		lassign $homeVectors($x,$y) hx hy
		# find closest neighbor
		lassign [closestNeighbor $x $y $hx $hy] xn yn dist
		# distance exceeded threshold
		if {$dist > $distThresh} {
		    # we have left the array
		    # mark all visited points as failure
		    # puts "$x $y failure"
		    set marker "$failure"
		    break
		}
		# draw line
		if {$flags & $pcaShowMovementLine} {
		    drawLine $canvas $x $y $xn $yn $movementColor 2
		}
		# move to next grid point
		set x $xn
		set y $yn
	    } else {
		# we have reached a point we have already visited
		# set marker accordingly
		set marker $state($x,$y)
		break
	    }
	    if {$flags & $pcaDemoMode} {
		update
		after $demoModeTimeDelay
	    }
	}
	# puts $visitedPoints
	# puts "marker: [lindex $stateChar $marker]"
	foreach pvis $visitedPoints {
	    lassign $pvis xvis yvis
	    set state($xvis,$yvis) $marker
	    if {$flags & $pcaMarkSuccessfulPoints} {
		# just mark successful points
		if {$marker == "$success"} {
		    drawCircle $canvas $xvis $yvis 0.15 black $markerColor
		}
	    }
	}
	if {$flags & $pcaDemoMode} {
	    update
	}
    }
    # print state array
    if {$flags & $pcaPrintStateArray} {
	for {set y $ymin} {$y <= $ymax} {incr y} {
	    for {set x $xmin} {$x <= $xmax} {incr x} {
		puts -nonewline [lindex $stateChar $state($x,$y)]
	    }
	    puts ""
	}
    }
    # compute RR
    # 27. Jun 12 (rm): changed computation: goal nodes are excluded
    set successCtr 0
    set totalCtr 0
    for {set x $xmin} {$x <= $xmax} {incr x} {
	for {set y $ymin} {$y <= $ymax} {incr y} {
	    if {$state($x,$y) != "$goal"} {
		incr totalCtr
	    } 
	    if {$state($x,$y) == "$success"} {
		incr successCtr
	    }
	}
    }
    set RR [expr double($successCtr) / $totalCtr]
    if {$flags & $pcaPlotWithGreyBoxes} {
	# grey values
	set black "\#000000"
	set lightgrey "\#C0C0C0"
	set darkgrey "\#808080"
	set grey($undecided) $black
	set grey($timeout) $lightgrey  
	set grey($failure) $lightgrey
	set grey($success) $darkgrey
	set grey($goal) $black
	# plot
	for {set x $xmin} {$x <= $xmax} {incr x} {
	    for {set y $ymin} {$y <= $ymax} {incr y} {
		# puts "$x $y $state($x,$y)"
		fillBox $canvas $x $y $grey($state($x,$y))
	    }
	}
    }
    finishPlot $canvas
    update
    return $RR
}

# ===========================================================================
# main
# ===========================================================================

# read input file
assert {[llength $argv] == 1}
set fn [lindex $argv 0]
puts "reading $fn"
set f [open $fn "r"]
assert {[gets $f ssPos] >= 0}
lassign $ssPos xss yss
puts "snapshot position ($xss, $yss)"
assert {[gets $f empty] == 0}
assert {[gets $f empty] == 0}
set xList [list]
set yList [list]
while {[gets $f homeVector] >= 0} {
    lassign $homeVector x y hx hy
    # need to normalize, deviations were too large after reading from file
    set hl [expr hypot($hx, $hy)]
    set hx [expr $hx / $hl]
    set hy [expr $hy / $hl]
    set homeVectors($x,$y) [list $hx $hy]
    lappend xList $x
    lappend yList $y
}
set xmin [tcl::mathfunc::min {*}$xList]
set xmax [tcl::mathfunc::max {*}$xList]
set ymin [tcl::mathfunc::min {*}$yList]
set ymax [tcl::mathfunc::max {*}$yList]
puts "grid xmin $xmin xmax $xmax ymin $ymin ymax $ymax"
close $f

# create canvas
createCanvas .c "Catchment Area"
wm iconify .
update
puts "plotting"
set RR [plotCatchmentArea .c $xss $yss $distThresh $goalMode $flags]
puts "return ratio $RR"
savePlot .c [file rootname $fn]_catch.ps
puts "done"
