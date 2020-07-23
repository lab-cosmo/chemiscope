import * as childProcess from 'child_process';
import * as fs from 'fs';
import * as path from 'path';
import * as tmp from 'tmp';

const tmpdir = tmp.dirSync({
    prefix: '.tmp',
    tmpdir: path.join(__dirname, '..'),
    unsafeCleanup: true,
});

childProcess.execSync(
    'git clone https://github.com/cosmo-epfl/chemiscope' +
        ' --depth=1 ' +
        ' --branch=gh-pages ' +
        tmpdir.name
);

for (const file of [
    'CSD-1000R.json.gz',
    'Arginine-Dipeptide.json.gz',
    'Qm7b.json.gz',
    'Azaphenacenes.json.gz',
    'Zeolites.json.gz',
]) {
    fs.renameSync(`${tmpdir.name}/${file}`, `./app/${file}`);
}

fs.mkdirSync('./app/structures/', { recursive: true });
for (let i = 0; i < 523; i++) {
    const file = `Azaphenacenes-${i}.json`;
    fs.renameSync(`${tmpdir.name}/structures/${file}`, `./app/structures/${file}`);
}

tmpdir.removeCallback();
