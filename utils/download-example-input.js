const fs = require('fs');
const url = require('url');
const https = require('https');

function download(file_url) {
    console.log(`downloading ${file_url}`);

    const basename = url.parse(file_url).pathname.split('/').pop();
    const file = fs.createWriteStream(`./app/${basename}`);

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
