import pytest
from nufebmgr.NufebProject import NufebProject

def test_initialization():
    project = NufebProject()
    assert project is not None

def test_error_on_assign_taxa_not_all_taxa_have_entries():
    prj = NufebProject()
    with pytest.raises(ValueError) as excinfo:
        prj.add_taxon_by_jsonfile("data/example_het_taxa.json")
        prj.layout_uniform(nbugs=20)
        prj.set_composition({'anammox_two_pathway': '33.333',
                             'denitrifier': '33.33',
                             'imperfect_denitrifier_no': '33.33'})
        prj.distribute_spatially_even()
        prj._assign_taxa()
    assert f'Setting composition with unknown taxon: anammox_two_pathway, denitrifier, imperfect_denitrifier_no' in str(excinfo.value)

    prj = NufebProject()
    try:
        prj.add_taxon_by_jsonfile("data/example_nitrogen_anaerobic_taxa.json")
        prj.layout_uniform(nbugs=20)
        prj.set_composition({'anammox_two_pathway': '33.333',
                             'denitrifier': '33.33',
                             'imperfect_denitrifier_no': '33.33'})
        prj.distribute_spatially_even()
        prj._assign_taxa()
    except Exception as e:
        pytest.fail(f'Unexpected exception {e}')

    prj = NufebProject()
    with pytest.raises(ValueError) as excinfo:
        prj.add_taxon_by_jsonfile("data/example_nitrogen_anaerobic_taxa_missing_amx.json")
        prj.layout_uniform(nbugs=20)
        prj.set_composition({'anammox_two_pathway': '33.333',
                             'denitrifier': '33.33',
                             'imperfect_denitrifier_no': '33.33'})
        prj.distribute_spatially_even()
        prj._assign_taxa()
    assert f'Setting composition with unknown taxon: anammox_two_pathway' in str(excinfo.value)

    prj = NufebProject()
    with pytest.raises(ValueError) as excinfo:
        prj.add_taxon_by_jsonfile("data/example_nitrogen_anaerobic_taxa_missing_amx_extra_het.json")
        prj.layout_uniform(nbugs=20)
        prj.set_composition({'anammox_two_pathway': '33.333',
                             'denitrifier': '33.33',
                             'imperfect_denitrifier_no': '33.33'})
        prj.distribute_spatially_even()
        prj._assign_taxa()
    assert f'Setting composition with unknown taxon: anammox_two_pathway' in str(excinfo.value)