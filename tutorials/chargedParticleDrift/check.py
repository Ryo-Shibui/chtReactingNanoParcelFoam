#!/usr/bin/env python3
"""Check a native charged-particle trajectory against its case gas properties."""
from pathlib import Path
import json,math,re
import numpy as np

ROOT = Path(__file__).resolve().parent

def list_payload(path):
    text = path.read_text().split('}', 1)[1]
    match = re.search(r'\b(\d+)\s*\((.*)\)', text, re.S)
    return int(match[1]), match[2]

def vectors(path):
    count, text = list_payload(path)
    array = np.asarray([
        list(map(float, value))
        for value in re.findall(
            r'\(\s*([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s*\)',
            text
        )
    ])
    assert len(array) == count
    return array

def scalars(path):
    count, text = list_payload(path)
    array = np.fromstring(text.replace('(', ' ').replace(')', ' '), sep=' ')
    assert len(array) == count
    return array

def scalar(text, key):
    return float(re.search(r'\b' + key + r'\s+([\deE.+-]+)\s*;', text)[1])

def dictionary(text, name):
    text = re.sub(r'//[^\n]*', '', text)
    start = re.search(r'\b' + name + r'\s*\{', text).end()
    index, depth = start, 1
    while depth:
        if text[index] == '{':
            depth += 1
        elif text[index] == '}':
            depth -= 1
        index += 1
    return text[start:index - 1]

def main():
    single = (ROOT / 'constant/foam.dat').exists()
    region = Path('.') if single else Path('fluid')
    constant = ROOT / 'constant' / region
    zero = ROOT / '0' / region
    nitrogen = dictionary((constant / 'foam.dat').read_text(), 'N2')
    cloud = (constant / 'reactingCloud1Properties').read_text()
    temperature = float(
        re.search(r'internalField\s+uniform\s+([^;]+);', (zero / 'T').read_text())[1]
    )
    mu = scalar(nitrogen, 'As') * temperature**1.5 / (
        temperature + scalar(nitrogen, 'Ts')
    )
    diameter = scalar(dictionary(cloud, 'fixedValueDistribution'), 'value')
    force = dictionary(cloud, 'nanoStokesDrag')
    meanFreePath = scalar(force, 'lambda' if re.search(r'\blambda\s', force) else 'meanFreePath')
    charge = scalar(cloud, 'parcelCharge')
    rhoP = scalar(cloud, 'rho0')
    electric = np.fromstring(
        re.search(r'internalField\s+uniform\s*\(([^)]+)\)', (zero / 'E').read_text())[1],
        sep=' '
    )
    kn = 2 * meanFreePath / diameter
    cc = 1 + kn * (1.257 + 0.4 * math.exp(-1.1 / kn))
    beta = 3 * math.pi * mu * diameter / cc
    mass = rhoP * math.pi * diameter**3 / 6
    tau = mass / beta
    terminal = charge * 1.602176634e-19 * electric / beta
    initial = vectors(constant / 'particlePositions')[0]
    rows = []
    for time in sorted(
        [p for p in ROOT.iterdir() if p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?', p.name)],
        key=lambda p: float(p.name)
    ):
        directory = time / region / 'lagrangian/reactingCloud1'
        if not (directory / 'positions').exists():
            continue
        t = float(time.name)
        position = vectors(directory / 'positions')[0]
        theoretical = initial + terminal * (t - tau * (1 - math.exp(-t / tau)))
        rows.append((t, *position, *theoretical))
    assert rows, 'No native parcel output was written'
    final = np.asarray(rows[-1])
    relative = np.linalg.norm(final[1:4] - final[4:7]) / np.linalg.norm(final[4:7] - initial)
    directory = ROOT / f'{final[0]:g}' / region / 'lagrangian/reactingCloud1'
    finalVelocity = vectors(directory / 'U')[0]
    velocityError = np.linalg.norm(finalVelocity - terminal) / np.linalg.norm(terminal)
    finalCharge = float(scalars(directory / 'chargeNumber')[0])
    passed = bool(relative < 0.002 and velocityError < 1e-9 and finalCharge == charge)
    output = ROOT / 'postProcessing/drift'
    output.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        output / 'trajectory.csv', np.asarray(rows), delimiter=',',
        header='time_s,x_m,y_m,z_m,theory_x_m,theory_y_m,theory_z_m', comments=''
    )
    result = dict(
        passed=passed, temperature_K=temperature, mu_N2_Pa_s=mu,
        terminalVelocity_m_s=terminal.tolist(), relaxationTime_s=tau,
        displacementRelativeError=float(relative), velocityRelativeError=float(velocityError),
        chargeNumber=finalCharge
    )
    (output / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if passed else 1)

if __name__ == '__main__':
    main()
