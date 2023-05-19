# vim: set expandtab ts=4 sw=4:

from chimerax.core.commands import (
    CmdDesc, FloatArg, Bounded, IntArg, BoolArg, FileNameArg
)
from chimerax.map import MapArg
from chimerax.atomic import ResiduesArg


def qscore(session, residues, to_volume=None, use_gui=True, reference_gaussian_sigma=0.6, points_per_shell=8,
            max_shell_radius=2.0, shell_radius_step=0.1, include_hydrogens=False, randomize_shell_points=True,
            log_details=False, output_file=None, map_resolution=3.0, window=0,
            show_res_in_map=None, show_atoms_in_map=None, dist_near=None):


    if show_res_in_map != None :
        print ( "showing res in map" )
        print ( show_res_in_map )
        from chimerax.atomic import Atoms
        resAtoms = Atoms()
        for r in residues :
            resAtoms = resAtoms + r.atoms
        from .other import showAtoms, MapFromId
        showAtoms ( session, None, resAtoms, dist_near, show_res_in_map, only=True )
        return

    if show_atoms_in_map != None :
        print ( "showing atoms in map" )
        print ( show_atoms_in_map )
        from chimerax.atomic import Atoms
        atoms = Atoms ( [] )
        for r in residues :
            atoms = atoms + r.atoms
        from .other import showAtoms, MapFromId
        showAtoms ( session, atoms, None, dist_near, show_atoms_in_map, only=True )
        return


    if to_volume is None:
        from chimerax.core.errors import UserError
        raise UserError("Must specify a map to compare the model to!")

    if use_gui :
        if not session.ui.is_gui:
            session.logger.warning('"qscore" command was called with "useGui True", but ChimeraX is in non-GUI mode. This argument has been ignored.')
        else:
            session.logger.warning('When the "qscore" command is called with "useGui True", the analysis will be run on the entire model.')
            from chimerax.core.commands import run
            from chimerax.qscore.tool import QScore_ToolUI
            run(session, 'ui tool show "Model-map Q-Score"')
            tool = session.tools.find_by_class(QScore_ToolUI)[0]
            mw = tool.tool_window.main_widget
            us = residues.unique_structures
            if not len(us):
                raise UserError('Selection contains no atoms!')
            if len(us) > 1:
                raise UserError('All residues should come from a single model!')
            mw.selected_model = us[0]
            mw.selected_volume = to_volume
            mw.reference_sigma = reference_gaussian_sigma
            mw.points_per_shell = points_per_shell
            mw.max_shell_radius = max_shell_radius
            mw.shell_radius_step = shell_radius_step
            mw.log_details = log_details
            residue_map, (query_atoms, atom_scores) = mw.recalc(log_details=log_details, output_file=output_file, echo_command=False)

    else :
        from .qscore import q_score
        residue_map, (query_atoms, atom_scores) = q_score(residues, to_volume, ref_sigma=reference_gaussian_sigma,
                                        points_per_shell=points_per_shell,
                                        max_rad=max_shell_radius,
                                        step=shell_radius_step,
                                        include_h=include_hydrogens,
                                        randomize_shell_points=randomize_shell_points,
                                        logger=session.logger,
                                        log_details=log_details,
                                        output_file=output_file,
                                        map_resolution=map_resolution,
                                        window=window
                                        )

    if not use_gui:
        # The GUI itself runs this command with use_gui=False, so if we don't do this the result gets printed twice
        session.logger.info(f'Overall mean Q-Score: {atom_scores.mean():.2f}')

    return residue_map, (query_atoms, atom_scores)



qscore_desc = CmdDesc(
    synopsis="Calculate the map-model Q-scores for a selection of residues",
    required=[
        ("residues", ResiduesArg),
    ],
    keyword=[
        ("to_volume", MapArg),
        ("use_gui", BoolArg),
        ("reference_gaussian_sigma", Bounded(FloatArg, min=0.1, max=2.0)),
        ("points_per_shell", Bounded(IntArg, min=2, max=32)),
        ("max_shell_radius", Bounded(FloatArg, min=0.5, max=2.5)),
        ("shell_radius_step", Bounded(FloatArg, min=0.025, max=0.5)),
        ("log_details", BoolArg),
        ("output_file", FileNameArg),
        ("map_resolution", Bounded(FloatArg, min=0.1, max=1000.0)),
        ("window", Bounded(IntArg, min=0, max=5)),
        ("show_res_in_map", MapArg),
        ("show_atoms_in_map", MapArg),
        ("dist_near", FloatArg),
    ]
    )




#

# toolshed uninstall chimerax-qscore; toolshed install /Users/greg/git/chimerax-qscore/dist/ChimeraX_QScore-1.0-cp39-cp39-macosx_10_13_universal2.whl
