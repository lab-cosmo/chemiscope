const fs = require('fs');
const url = require('url');
const https = require('https');

function download(file_url) {
    console.log(`downloading ${file_url}`);

    const pathname = url.parse(file_url).pathname;
    const file = fs.createWriteStream(`./app/${pathname}`);

    const req = https.request(file_url, res => {
        res.on('data', data => file.write(data))
           .on('end', () => file.end());
    });
    req.end();
}

download('https://chemiscope.org/CSD-500.json.gz');
download('https://chemiscope.org/Arginine-Dipeptide.json.gz');
download('https://chemiscope.org/Qm7b.json.gz');
download('https://chemiscope.org/Azaphenacenes.json.gz');

fs.mkdirSync('./app/structures/', {recursive: true});
for (let i=0; i<523; i++) {
    // separate requests to not get banned by DoS protection
    const delay = (i % 100) * 100;
    setTimeout(() => download(`https://chemiscope.org/structures/Azaphenacenes-${i}.json`), delay);
}
