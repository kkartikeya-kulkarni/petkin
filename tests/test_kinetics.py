import numpy as np
import pytest
from scipy.integrate import solve_ivp
from petkin import Frames, InputFunction, predict, fit_compartment, logan


def fixture():
    t=np.linspace(0,90,9001)
    inp=InputFunction(t, 8*t*np.exp(-t/.7)+1.5*np.exp(-t/25))
    edges=np.r_[np.arange(0,5,.25), np.arange(5,20,1), np.arange(20,91,2)]
    return Frames(edges[:-1],edges[1:]),inp


def independent(frames,inp,p):
    """Adaptive ODE + integrated observation state; independent of convolution code."""
    K1,k2,k3,k4=p
    def rhs(t,y):
        cp=np.interp(t,inp.time,inp.plasma)
        return [K1*cp-(k2+k3)*y[0]+k4*y[1],k3*y[0]-k4*y[1],y[0]+y[1]]
    sol=solve_ivp(rhs,(0,frames.end[-1]),[0,0,0],rtol=1e-9,atol=1e-11,max_step=.03,dense_output=True)
    return (sol.sol(frames.end)[2]-sol.sol(frames.start)[2])/frames.duration


@pytest.mark.parametrize('model,p', [('1t',[.17,.22]),('2t',[.17,.22,.09,.065]),('2t_irreversible',[.17,.22,.09])])
def test_independent_recovery(model,p):
    frames,inp=fixture()
    q=list(p)+([0,0] if model=='1t' else [0] if model=='2t_irreversible' else [])
    y=independent(frames,inp,q)
    np.testing.assert_allclose(predict(frames,inp,p,model,dt=.01),y,rtol=.003,atol=1e-4)
    fit=fit_compartment(frames,inp,y,model,dt=.01)
    assert fit['success']
    np.testing.assert_allclose(list(fit['parameters'].values()),p,rtol=.01)


def test_logan_1t():
    f,inp=fixture(); y=independent(f,inp,[.2,.15,0,0])
    out=logan(f,y,inp,30)
    assert out['VT']==pytest.approx(.2/.15,rel=.015)
    assert out['r_squared']>.999


def test_blood_and_convergence():
    f,inp=fixture(); inp=InputFunction(inp.time,inp.plasma,inp.plasma*.8)
    p=[.2,.15]
    coarse=predict(f,inp,p,'1t',vb=.06,dt=.025)
    fine=predict(f,inp,p,'1t',vb=.06,dt=.005)
    np.testing.assert_allclose(coarse,fine,rtol=.002,atol=1e-4)
    assert logan(f,fine,inp,30,vb=.06)['VT']==pytest.approx(.2/.15,rel=.015)


def test_invalid_inputs():
    with pytest.raises(ValueError): Frames([0,1],[2,3])
    with pytest.raises(ValueError): InputFunction([1,2],[1,2])
    f,inp=fixture()
    with pytest.raises(ValueError): predict(f,inp,[.1,.2],'1t',vb=.1)
    with pytest.raises(ValueError): fit_compartment(f,inp,np.ones(len(f.start)),sigma=np.zeros(len(f.start)))
    with pytest.raises(ValueError): logan(Frames([1,2,3],[2,3,4]),[1,2,3],inp,1)
    with pytest.raises(ValueError): logan(f,np.zeros(len(f.start)),inp,30)


def test_frame_boundary_roundoff():
    starts=np.array([126,138,150])/60
    duration=np.array([12,12,30])/60
    Frames(starts,starts+duration)
