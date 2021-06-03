import * as childProcess from 'child_process';
import * as fs from 'fs';
import * as path from 'path';
import * as tmp from 'tmp';

const ALL_EXAMPLES = [
    'CSD-1000R.json.gz',
    'Arginine-Dipeptide.json.gz',
    'Qm7b.json.gz',
    'Azaphenacenes.json.gz',
    'Zeolites.json.gz',
];

let needsUpdate = false;
for (const file of ALL_EXAMPLES) {
    if (!fs.existsSync(`./app/examples/${file}`)) {
        needsUpdate = true;
    }
}

if (!needsUpdate) {
    // eslint-disable-next-line no-console
    console.log(
        'Example input files already exists. Run `rm -rf app/examples/*.json.gz` if you want to update them'
    );
    process.exit();
}

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

for (const file of ALL_EXAMPLES) {
    fs.renameSync(`${tmpdir.name}/examples/${file}`, `./app/examples/${file}`);
}

fs.mkdirSync('./app/examples/structures/', { recursive: true });
for (let i = 0; i < 523; i++) {
    const file = `Azaphenacenes-${i}.json`;
    fs.renameSync(
        `${tmpdir.name}/examples/structures/${file}`,
        `./app/examples/structures/${file}`
    );
}

tmpdir.removeCallback();
