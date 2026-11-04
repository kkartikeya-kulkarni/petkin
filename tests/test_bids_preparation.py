import hashlib,json
import numpy as np
import pytest
from petkin.bids import input_from_bids,prepare_bids_input
from petkin.petprep import import_petprep_tacs
from petkin.snapshot import fetch_snapshot


def records(tmp_path):
    pet=dict(FrameTimesStart=[0,30],FrameDuration=[30,30],TimeZero='12:00:00',ScanStart=0,InjectionStart=0,Units='Bq/mL',ImageDecayCorrected=True,ImageDecayCorrectionTime=0,FrameReferenceTime=[9000,9999])
    blood=dict(time={'Units':'s'},plasma_radioactivity={'Units':'kBq/mL'},whole_blood_radioactivity={'Units':'Bq/mL'},WholeBloodAvail=True,MetaboliteAvail=True,DispersionCorrected=False)
    (tmp_path/'pet.json').write_text(json.dumps(pet));(tmp_path/'blood.json').write_text(json.dumps(blood))
    (tmp_path/'blood.tsv').write_text('time\tplasma_radioactivity\twhole_blood_radioactivity\tmetabolite_parent_fraction\n10\t4\t3000\t0.5\n40\t2\t1500\t0.25\n')
    return pet,blood


def test_total_plasma_corrected_once_whole_blood_unchanged_and_timezero_preserved(tmp_path):
    records(tmp_path);qc=prepare_bids_input(tmp_path/'blood.tsv',tmp_path/'blood.json',tmp_path/'pet.json',tmp_path/'out',plasma_kind='total',initial_policy='zero',terminal_policy='hold')
    result=np.genfromtxt(tmp_path/'out/input.csv',names=True,delimiter=',')
    np.testing.assert_allclose(result['time'],[0,10/60,40/60,1]);np.testing.assert_allclose(result['plasma'],[0,2000,500,500]);np.testing.assert_allclose(result['blood'],[0,3000,1500,1500])
    assert qc['time_origin']=='12:00:00' and qc['assumptions']==['initial_zero','terminal_hold']
    assert not qc['PET_blood_cross_calibration_verified'] and not qc['blood_decay_reference_verified']
    assert qc['dispersion_corrected'] is False


def test_processed_parent_plasma_is_not_corrected_twice(tmp_path):
    _,blood=records(tmp_path)
    inp,qc=input_from_bids(tmp_path/'blood.tsv',blood,plasma_mode='processed',initial_policy='zero')
    assert inp.plasma[1]==4000 and inp.blood[1]==3000

@pytest.mark.parametrize('fraction',['-0.1','1.1','nan'])
def test_invalid_parent_fraction_rejected(tmp_path,fraction):
    _,blood=records(tmp_path);p=tmp_path/'blood.tsv';p.write_text(p.read_text().replace('0.5',fraction))
    with pytest.raises(ValueError,match='Parent fractions'):input_from_bids(p,blood,plasma_mode='total',initial_policy='zero')


def test_missing_parent_fraction_never_guessed(tmp_path):
    _,blood=records(tmp_path);p=tmp_path/'blood.tsv';p.write_text('time\tplasma_radioactivity\twhole_blood_radioactivity\n10\t4\t3000\n40\t2\t1500\n')
    with pytest.raises(ValueError,match='parent fractions'):input_from_bids(p,blood,plasma_mode='total',initial_policy='zero')


def test_petprep_import_verifies_seconds_and_keeps_missing_region_values_visible(tmp_path):
    records(tmp_path);p=tmp_path/'tacs.tsv';p.write_text('frame_start\tframe_end\tCortex left\tThalamus\n0\t30\t10\tn/a\n30\t60\t15\t8\n')
    result=import_petprep_tacs(p,tmp_path/'pet.json',tmp_path/'regions')
    assert result['regions'][0]['name']=='Cortex left' and result['regions'][1]['nonfinite_frames']==[0]
    tac=np.genfromtxt(tmp_path/'regions/region_0000.csv',delimiter=',',names=True);np.testing.assert_allclose(tac['end'],[.5,1]);np.testing.assert_allclose(tac['activity'],[10,15])
    p.write_text(p.read_text().replace('30','0.5').replace('60','1'))
    with pytest.raises(ValueError,match='timing'):import_petprep_tacs(p,tmp_path/'pet.json',tmp_path/'bad')


def test_petprep_duplicate_regions_and_mismatched_frame_counts_rejected(tmp_path):
    records(tmp_path);p=tmp_path/'tacs.tsv';p.write_text('frame_start\tframe_end\tCortex\tCortex\n0\t30\t1\t2\n30\t60\t3\t4\n')
    with pytest.raises(ValueError,match='distinct'):import_petprep_tacs(p,tmp_path/'pet.json',tmp_path/'out')


def test_snapshot_fetch_keeps_pinned_metadata_and_checks_renamed_cached_image(tmp_path,monkeypatch):
    from petkin import snapshot
    metadata=b'{"TimeZero":"12:00:00"}';image=b'known image bytes';cache=tmp_path/'cache';cache.mkdir();(cache/'sub-OLD_pet.nii').write_bytes(image)
    commit='a'*40;gitsha=hashlib.sha1(b'blob '+str(len(metadata)).encode()+b'\0'+metadata).hexdigest()
    files=[dict(path='sub-NEW/pet.json',kind='git',git_blob_sha=gitsha),dict(path='sub-NEW/pet.nii',kind='annex',bytes=len(image),sha256=hashlib.sha256(image).hexdigest())]
    lock=tmp_path/'lock.json';lock.write_text(json.dumps(dict(dataset='ds004868',commit=commit,requested_revision='1.0.4',files=files)))
    urls=[]
    def read(url):urls.append(url);return metadata
    monkeypatch.setattr(snapshot,'_read',read)
    result=fetch_snapshot(lock,tmp_path/'out',cache)
    assert result['all_files_verified'] and (tmp_path/'out/sub-NEW/pet.nii').read_bytes()==image
    assert len(urls)==1 and commit in urls[0] and 's3.amazonaws' not in urls[0]
    monkeypatch.setattr(snapshot,'_read',lambda url:b'wrong metadata')
    with pytest.raises(ValueError,match='Git blob'):fetch_snapshot(lock,tmp_path/'corrupt',cache)
    assert not (tmp_path/'corrupt/sub-NEW/pet.json').exists()


def test_snapshot_rejects_traversal_paths_without_writing(tmp_path):
    lock=tmp_path/'lock.json';lock.write_text(json.dumps(dict(dataset='ds004868',commit='a'*40,requested_revision='1.0.4',files=[dict(path='../escape',kind='git',git_blob_sha='a'*40)])))
    with pytest.raises(ValueError,match='Unsafe'):fetch_snapshot(lock,tmp_path/'out')
    assert not (tmp_path/'out').exists()


def test_injection_before_origin_cannot_be_silently_treated_as_zero_initial_tissue(tmp_path):
    pet,_=records(tmp_path);pet['InjectionStart']=-12;(tmp_path/'pet.json').write_text(json.dumps(pet))
    with pytest.raises(ValueError,match='Injection predates TimeZero'):
        prepare_bids_input(tmp_path/'blood.tsv',tmp_path/'blood.json',tmp_path/'pet.json',tmp_path/'out',plasma_kind='total',initial_policy='zero',terminal_policy='hold')
    assert not (tmp_path/'out').exists()


def test_nonzero_pet_decay_reference_rejected_without_invented_rescaling(tmp_path):
    pet,_=records(tmp_path);pet['ImageDecayCorrectionTime']=20;(tmp_path/'pet.json').write_text(json.dumps(pet))
    with pytest.raises(ValueError,match='decay-corrected to TimeZero'):
        prepare_bids_input(tmp_path/'blood.tsv',tmp_path/'blood.json',tmp_path/'pet.json',tmp_path/'out',plasma_kind='total',initial_policy='zero',terminal_policy='hold')


def test_different_scan_inputs_rejected_before_analysis(tmp_path):
    from petkin.workflow import analyze_regions
    records(tmp_path)
    prepare_bids_input(tmp_path/'blood.tsv',tmp_path/'blood.json',tmp_path/'pet.json',tmp_path/'input',plasma_kind='total',initial_policy='zero',terminal_policy='hold')
    p=tmp_path/'tacs.tsv';p.write_text('frame_start\tframe_end\tCortex\n0\t30\t10\n30\t60\t15\n')
    meta=json.loads((tmp_path/'pet.json').read_text());meta['TracerName']='another acquisition';(tmp_path/'other.json').write_text(json.dumps(meta))
    import_petprep_tacs(p,tmp_path/'other.json',tmp_path/'regions')
    with pytest.raises(ValueError,match='different PET metadata'):analyze_regions(tmp_path/'regions/regions.json',tmp_path/'input/input.csv',tmp_path/'out')
    assert not (tmp_path/'out').exists()


def test_snapshot_lock_retains_reference_transforms_without_guessing_file_extensions(tmp_path,monkeypatch):
    from petkin import snapshot
    commit='b'*40
    tree={'tree':[{'path':'sub-X/pet/reference_transform.dat','type':'blob','mode':'100644','sha':'c'*40,'size':20}]}
    def read(url):
        if '/commits/' in url:return json.dumps({'sha':commit}).encode()
        return json.dumps(tree).encode()
    monkeypatch.setattr(snapshot,'_read',read)
    result=snapshot.lock_snapshot('ds004868','1.0.4',tmp_path/'lock.json',participant='X')
    assert result['commit']==commit and result['files'][0]['path'].endswith('reference_transform.dat')


def test_legacy_dicom_units_normalized_without_rewriting_pinned_metadata(tmp_path):
    from petkin.bids import validate_pet_metadata
    pet,_=records(tmp_path);pet['Units']='BQML';(tmp_path/'pet.json').write_text(json.dumps(pet));before=(tmp_path/'pet.json').read_bytes()
    checked=validate_pet_metadata(pet)
    assert checked['units']=='Bq/mL' and checked['source_units']=='BQML' and checked['unit_alias_applied']
    prepare_bids_input(tmp_path/'blood.tsv',tmp_path/'blood.json',tmp_path/'pet.json',tmp_path/'out',plasma_kind='total',initial_policy='zero',terminal_policy='hold')
    assert (tmp_path/'pet.json').read_bytes()==before
