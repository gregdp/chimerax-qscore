

##
##@Author: Greg Pintilie
##@Date:   14-Apr-2023
##@Email:  gregdp@gmail.com
##@Last modified by:   gregdp
##@Last modified time: 14-Apr-2023
##@License: Free for non-commercial use (see license.pdf)
##@Copyright: Copyright 2023 Grigore Pintilie
##



import chimerax

def listRes ( session, resNum=None ) :

    for m in session.models :
        if m.visible == True and type ( m ) == chimerax.atomic.structure.AtomicStructure :
            print ( m.name )

            for ri, res in enumerate ( m.residues ) :
                if resNum != None :
                    if resNum == ri+1 :
                        print ( "%s - %d - %d - %s" % (res.chain_id, ri+1, res.number, res.name) )
                else :
                    print ( "%s - %d - %d - %s" % (res.chain_id, ri+1, res.number, res.name) )


def ModId ( mol ) :
    mid = "%d" % mol.id[0]
    for i in range ( 1, len(mol.id) ) :
        mid += ".%d" % mol.id[i]
    return mid

def showLigands ( session ) :

    from chimerax.atomic import Sequence

    for m in session.models :
        # print ( "-- %s --" % m.name, type(m) )
        if m.visible and type(m) == chimerax.atomic.structure.AtomicStructure :
            modId = ModId(m)
            for res in m.residues :
                isProt = Sequence.protein3to1(res.name) != "X"
                isNA = Sequence.nucleic3to1(res.name) != "X"
                if not isProt and not isNA :
                    ss = "#%s/%s:%s" % (modId, res.chain_id, res.number)
                    print ( "Chain %s:%d %s/%s -- %s" % (res.chain_id, res.number, res.name, Sequence.protein3to1(res.name), ss) )



def showAtoms ( session, atoms, resAtoms, distNear, inMap, only=True ) :

    if atoms == None and resAtoms == None :
        print ( " - nothing to show?" )
        return

    #ress = session.selection.items ("residues")[0] # chimerax.atomic.molarray.Residues
    #print ( " - %d selected residues" % len(ress) )

    if atoms != None and len(atoms) == 0 :
        raise Exception ( "no atoms specified by atoms keyword - check model id?" )

    if resAtoms != None and len (resAtoms) == 0 :
        raise Exception ( "no atoms specified by res keyword - check model id?" )

    mol = atoms[0].structure if atoms != None else resAtoms[0].structure

    from chimerax.atomic import Atoms
    atoms = atoms if atoms != None else Atoms ( [] )
    resAtoms = resAtoms if resAtoms != None else Atoms ( [] )

    resm = {}
    for at in resAtoms :
        resm[at.residue] = 1

    print ( "Showing %d atoms, %d residues" % (len(atoms), len(resm)) )

    #if len(atoms) == 0 and len(resAtoms) == 0 :
    #    print ( "No atoms/residues specified" )


    amap = {}
    zatoms = {}
    for at in atoms :
        at.display = True
        at.draw_mode = at.STICK_STYLE
        amap[at] = 1
        zatoms[at] = 1

    for res in resm.keys() :
        if res.polymer_type == res.PT_PROTEIN or res.polymer_type == res.PT_NUCLEIC :
            res.ribbon_display = True
        else :
            for at in res.atoms :
                amap[at] = 1
                at.display = True
                at.draw_mode = at.STICK_STYLE

        for at in res.atoms :
            if at.element.number != 1 :
                zatoms[at] = 1

    if only :
        for at in mol.atoms :
            if not at in amap :
                at.display = False
                if not at.residue in resm :
                    at.residue.ribbon_display = False

        for res in mol.residues :
            if not res in resm :
                res.ribbon_display = False

    if distNear != None :
        print ( "Adding points within %f" % distNear )
        from chimerax.geometry import find_close_points, find_closest_points

        resm2 = {}
        for at in atoms :
            _, nearby_i = find_close_points( atoms.coords, mol.atoms.coords, distNear )
            for ai in nearby_i :
                resm2 [ mol.atoms[ai].residue ] = 1

        for res in resm2.keys () :
            for at in res.atoms :
                at.display = True
                at.draw_mode = at.STICK_STYLE
                if at.element.number != 1 :
                    zatoms[at] = 1
            if res.polymer_type == res.PT_PROTEIN or res.polymer_type == res.PT_NUCLEIC :
                res.ribbon_display = False


        resm2 = {}
        for res in resm.keys() :
            _, nearby_i = find_close_points( res.atoms.coords, mol.atoms.coords, distNear )
            for ai in nearby_i :
                resm2 [ mol.atoms[ai].residue ] = 1

        for res in resm2.keys () :
            for at in res.atoms :
                at.display = True
                at.draw_mode = at.STICK_STYLE
                if at.element.number != 1 :
                    zatoms[at] = 1
            if res.polymer_type == res.PT_PROTEIN or res.polymer_type == res.PT_NUCLEIC :
                res.ribbon_display = True

    closem = []
    for m in session.models :
        if m.name == "Q_scores zone" :
            closem.append ( m )
    session.models.close ( closem )

    if inMap != None :
        print ( "Adding zone map: %s" % inMap.name  )
        ZoneMapWithAtoms ( session, Atoms(zatoms.keys()), inMap, 2.0, "Q_scores zone" )
        inMap.display = False




def ZoneMapWithAtoms ( session, atoms, dmap, atomRad, nname, showMesh = False, alpha=0.4 ) :

    from numpy import min, max, zeros, float32, single, copy

    points1 = atoms.scene_coords
    points1 = dmap.position.inverse() * points1
    points0 = copy ( points1 )
    points1 = dmap.data.xyz_to_ijk_transform * points1

    bound = 5
    li,lj,lk = min ( points1, axis=0 ) - (bound, bound, bound)
    hi,hj,hk = max ( points1, axis=0 ) + (bound, bound, bound)

    n1 = hi - li + 1
    n2 = hj - lj + 1
    n3 = hk - lk + 1

    #print ( " - bounds: %d %d %d --> %d %d %d --> %d %d %d" % ( li,lj,lk, hi,hj,hk, n1,n2,n3 ) )
    nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )
    nn1 = int ( round (dmap.data.step[0] * float(n1) / nstep[0]) )
    nn2 = int ( round (dmap.data.step[1] * float(n2) / nstep[1]) )
    nn3 = int ( round (dmap.data.step[2] * float(n3) / nstep[2]) )

    O = dmap.data.origin
    #print " - %s origin:" % dmap.name, O
    nO = ( O[0] + float(li) * dmap.data.step[0],
           O[1] + float(lj) * dmap.data.step[1],
           O[2] + float(lk) * dmap.data.step[2] )
    #print ( " - new map origin:", nO )

    nmat = zeros ( (nn1,nn2,nn3), float32 )
    #ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles )
    from chimerax.map_data import ArrayGridData
    ndata = ArrayGridData ( nmat, nO, nstep )

    from chimerax.map_data import grid_indices
    npoints = grid_indices ( (nn1, nn2, nn3), single)  # i,j,k indices
    npoints = ndata.ijk_to_xyz_transform * npoints

    dvals = dmap.interpolated_values ( npoints ) # in local volume coordinates
    nmat = dvals.reshape( (nn3,nn2,nn1) )
    ndata = ArrayGridData(nmat, nO, nstep)

    from chimerax.map import volume_from_grid_data
    #v = volume_from_grid_data(ndata, session, style = 'surface', open_model = False, show_dialog = True)
    #levels, colors = dmap.initial_surface_levels()
    #v.set_parameters(surface_levels = levels, surface_colors = colors)
    #session.models.add([v])
    #nv.openState.xform = dmap.openState.xform

    from chimerax.map_data import zone_masked_grid_data
    mdata = zone_masked_grid_data(ndata, points0, atomRad)
    v = volume_from_grid_data(mdata, session, style = 'surface', open_model = False, show_dialog = True)
    levels, colors = dmap.initial_surface_levels()

    ct = []
    if len(dmap.surfaces) > 0 :
        levels = ( [dmap.surfaces[0].level] )
        C = dmap.surfaces[0].color
        ct = [ [C[0]/255.0, C[1]/255.0, C[2]/255.0, alpha] ]
    else :
        for c in colors :
            ct.append ( [c[0], c[1], c[2], alpha] )

    v.set_parameters(surface_levels = levels, surface_colors = ct)
    v.name = nname
    session.models.add([v])

    return v

    ro = VolumeViewer.volume.Rendering_Options()
    ro.smoothing_factor = .3
    ro.smoothing_iterations = 2
    ro.surface_smoothing = False
    ro.square_mesh = True
    ro.line_thickness = 2
    nv.update_surface ( False, ro )
    setro (ro)
    for sp in nv.surfacePieces :
        v, t = sp.geometry
        if len(v) == 8 and len(t) == 12 :
            sp.display = False
        else :
            if showMesh :
                sp.color = (.5, .5, .5, 1.0)
                sp.displayStyle = sp.Mesh
            else :
                sp.color = (0.7, 0.7, 0.7, alpha)

    return nv



def MapFromId ( session, inMap ) :
    dmap = None
    print ( "looking for" )
    print ( inMap )
    if inMap != None :
        #print ( inMap, len(inMap) )
        for mod in session.models :
            if mod.id == inMap :
                print ( " - found map: %s" % mod.name )
                print ( type(mod) )
                if type(mod) == chimerax.map.volume.Volume :
                    return mod
        raise Exception ( "inMap parameter did not find a map - check model id?" )



# dav res resNum 646
# swap /E:96 CYS criteria d density #3
# swap /E:96 CYS criteria c











#
