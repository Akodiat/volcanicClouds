import * as THREE from "three";
import {MapControls} from "./libs/threeAddons/MapControls.js";
import {Lut} from "./libs/threeAddons/Lut.js";
import ThreeGeo from "./libs/three-geo-esm.js";
import {tomoInverse} from "./src/tomoInverse.js";
import {locateVolcano} from "./src/locateVolcano.js";
import {drawParticles} from "./libs/draw.js";
import {TomographicPlaneGeometry} from "./src/tomographicPlaneGeometry.js";
import {GLTFLoader} from "./libs/threeAddons/GLTFLoader.js";
import {GLTFExporter} from "./libs/threeAddons/GLTFExporter.js";

let camera, scene, renderer, controls;
let frames = [];
let currentFrame;
init();
render();

// Initialise scene
function init() {
    // Setup renderer
    renderer = new THREE.WebGLRenderer({
        alpha: true
    });
    renderer.setPixelRatio(window.devicePixelRatio);
    renderer.setSize(window.innerWidth, window.innerHeight);
    renderer.toneMapping = THREE.ACESFilmicToneMapping;
    const container = document.getElementById("container");
    container.appendChild(renderer.domElement);

    // Setup scene and camera
    scene = new THREE.Scene();
    camera = new THREE.PerspectiveCamera(55, window.innerWidth / window.innerHeight, 0.01, 1e4);

    // Setup lights
    const pointLight = new THREE.PointLight(0xff0000, 100);
    pointLight.position.set(1, 1, 1);
    scene.add(pointLight);

    // Add x-y-z axis indicator
    const axesHelper = new THREE.AxesHelper(5);
    scene.add(axesHelper);

    // And camera controls
    controls = new MapControls(camera, renderer.domElement);
    controls.maxPolarAngle = Math.PI/2;
    controls.zoomToCursor = true;
    controls.addEventListener("change", render);

    // Load data when file is uploaded
    const fileInput = document.getElementById("fileInput");
    const loadFromFiles = async () => {
        const parseDate = (d, t) => new Date(
            `${d.slice(0,4)}-${d.slice(4,6)}-${d.slice(6,8)}T${t.slice(0,2)}:${t.slice(2,4)}`
        );
        const data = [];
        let alreadyProcessedData = [];
        for (const file of fileInput.files) {
            const text = await file.text();
            const [filename, suffix] = file.name.split(".");
            if (suffix === "csv") {
                const frame = {points: []};

                // Parse datetime
                const [_0, _1, day1, time1, _2, day2, time2] = filename.split("_");
                const date1 = parseDate(day1, time1);
                const date2 = parseDate(day2, time2);
                // Calc average time
                frame.time = new Date((
                    date1.getTime() +
                    date2.getTime()
                ) / 2);
                // Parse points
                for (let line of text.split("\n")) {
                    line = line.trim("\r");
                    if (line === "") {
                        continue;
                    }
                    const values = line.split(",").map(v=>parseFloat(v));
                    if (values.length == 2) {
                        [frame.size1, frame.size2] = values;
                    } else {
                        const [lonPutm, latPutm, altP, Concentration] = values;
                        frame.points.push({lonPutm, latPutm, altP, Concentration});
                    }
                }
                alreadyProcessedData.push(frame);
            } else {
                const scans = parseScans(text);
                data.push(scans);
            }
        }
        document.getElementById("fileUploadContainer").style.display = "none";
        onDataLoaded(data, alreadyProcessedData);
    };

    fileInput.onchange = loadFromFiles;

    window.addEventListener("keydown", (event) => {
        switch (event.code) {
            case "Enter":
                if (fileInput.files.length > 0) {
                    loadFromFiles();
                }
                break;
        }
    });

    // Update camera aspect ratio on window resize
    window.addEventListener("resize", onWindowResize);

    render();
}

function parseScans(text) {
    // Match spectral header and data
    const rInfo = /<scaninformation>(?<info>([\s\S])*?)<\/scaninformation>/gm;
    const scans = [];
    for (const match of text.matchAll(rInfo)) {
        const scanInfo = {};
        const info = match.groups.info.split("\n");
        info.forEach(d=>{
            d = d.trim("\r");
            if (d !== "") {
                const [k, v] = d.split("=");
                scanInfo[k] = isNaN(v) ? v : Number(v);
            }
        });
        scans.push({
            scanInfo: scanInfo
        });
    }

    // Match spectral header and data
    const rData = /#(?<header>[\s\S]+?)<spectraldata>(?<data>([\s\S])*?)<\/spectraldata>/gm;
    //const scans = [];
    let i = 0;
    for (const match of text.matchAll(rData)) {
        const spectralData = [];
        // Header names are not consistent across different stations
        // So let's assume that the order is at least the same
        //const header = match.groups.header.split("\t");
        const header = ["scanangle", "starttime", "stoptime", "name", "specsaturation", "fitsaturation", "counts_ms", "delta", "chisquare", "exposuretime", "numspec", "column_SO2", "columnerror_SO2", "shift_SO2", "shifterror_SO2", "squeeze_SO2", "squeezeerror_SO2", "column_O3", "columnerror_O3", "shift_O3", "shifterror_O3", "squeeze_O3", "squeezeerror_O3", "column_RING", "columnerror_RING", "shift_RING", "shifterror_RING", "squeeze_RING", "squeezeerror_RING", "isgoodpoint", "offset", "flag"];
        const data = match.groups.data.split("\n");
        data.forEach(d=>{
            d = d.trim("\r");
            if (d !== "") {
                const linedata = {};
                d.split("\t").forEach((v, i) => {
                    linedata[header[i].trim()] = isNaN(v) ? v : Number(v);
                });
                spectralData.push(linedata);
            }
        });
        scans[i].spectralData = spectralData;
        i++;
    }
    return scans;
}

function onDataLoaded(data, processedData) {
    console.log(data);

    let tgeo = new ThreeGeo();

    const [nameVol, latVol, lonVol, altVol] = locateVolcano(data);
    const summitLatLng = new THREE.Vector2(latVol, lonVol);
    const radius = 6.0;

    const loader = new GLTFLoader().setPath('resources/terrainMeshes/');

    const filename = `${nameVol}.glb`
    loader.load(filename, gltf => {
        const model = gltf.scene;
        scene.add(model);
        render();
    }, undefined, ()=>{
        // On error (file not found)
        const tokenMapbox = prompt(`The terrain for the volcano ${nameVol} is not saved. Input a mapbox token to download. To avoid this in the future, save the downloaded file to ./resources/terrainMeshes/`);
        tgeo = new ThreeGeo({
            tokenMapbox: tokenMapbox,
        });
        tgeo.getTerrainRgb(
            summitLatLng.toArray(),  // [lat, lng]
            radius,            // radius of bounding circle (km)
            13                 // zoom resolution
        ).then(terrain => {
            terrain.rotation.x = - Math.PI/2;
            scene.add(terrain);
            render();

            const gltfExporter = new GLTFExporter();
            gltfExporter.parse(
                terrain,
                function (result) {
                    saveArrayBuffer(result, filename);
                },
                error => console.log('An error happened during parsing', error),
                {binary: true}
            );
        });
    } );


    const {proj, unitsPerMeter} = tgeo.getProjection(summitLatLng.toArray(), radius);

    const toSceneCoords = (latLng, altitude) => {
        const pos2D = new THREE.Vector2(...proj(latLng));
        return new THREE.Vector3(pos2D.x, altitude * unitsPerMeter, -pos2D.y);
    };

    const summitPos = toSceneCoords(summitLatLng, altVol, proj);

    controls.minDistance = unitsPerMeter;
    controls.target.copy(summitPos);
    controls.update();

    const instPos = [];
    for (const instrumentData of data) {
        const scanInfo = instrumentData[0].scanInfo; // Use first datapoint
        const instrumentLatLng = new THREE.Vector2(
            scanInfo.lat,
            scanInfo.long
        );
        const instrumentPos = toSceneCoords(instrumentLatLng, scanInfo.alt, proj);
        instPos.push(instrumentPos);

        // Add a cube
        const instrumentGeometry = new THREE.BoxGeometry(1.5, 1, 3);
        const instrumentMaterial = new THREE.MeshStandardMaterial({color: 0x00ff00});
        const cube = new THREE.Mesh(instrumentGeometry, instrumentMaterial);
        cube.scale.multiplyScalar(100 * unitsPerMeter);
        cube.position.copy(instrumentPos);
        cube.lookAt(summitPos);
        scene.add(cube);
    }

    let positions = [];
    if (processedData.length === 0) {
        const deg2utm = (lat, long) => {
            const [x,y] = proj([lat, long]);
            return [x/unitsPerMeter, y/unitsPerMeter];
        }
        processedData = tomoInverse(data, deg2utm);
        for (const frame of processedData) {
            positions.push(frame.points.map(d=>new THREE.Vector3(
                    d.latPutm * unitsPerMeter,
                    d.altP * unitsPerMeter,
                    d.lonPutm * unitsPerMeter,
                    proj
            )));
        }

    } else {
        for (const frame of processedData) {
            positions.push(frame.points.map(d=>toSceneCoords(
                new THREE.Vector2(
                    d.latPutm, d.lonPutm
                ), d.altP, proj
            )));
        }
    }
    let dir;
    let t = 0;
    for (const frame of processedData) {
        const concentrations = frame.points.map(d=>d.Concentration);
        const lut = new Lut("ylOrRd", 512);
        lut.minV = Math.min(...concentrations);
        lut.maxV = Math.max(...concentrations);
        const colors = concentrations.map(c=>{
            const color = lut.getColor(c);
            return color;
        });
        const ps = positions[t];

        // Particles

        const particles = drawParticles(ps, colors, 0.005);

        // Tomographic plane

        const texture = new THREE.CanvasTexture(
            generateTexture(concentrations, frame.size1, frame.size2)
        );
        texture.wrapS = THREE.ClampToEdgeWrapping;
        texture.wrapT = THREE.ClampToEdgeWrapping;
        texture.colorSpace = THREE.SRGBColorSpace;

        const material = new THREE.MeshBasicMaterial({
            side: THREE.DoubleSide,
            map: texture,
            transparent: true
        });

        const line = instPos[0].clone().sub(
            instPos[1]
        );
        dir = line.clone().cross(new THREE.Object3D().up);
        if ((instPos[0].lengthSq() > instPos[0].clone().add(dir).lengthSq())) {
            // Make sure dir points away from the volcano
            // The volcano is at the origin,
            // so length gets the distance to it
            dir.negate();
        }

        const planeGeometry = new TomographicPlaneGeometry(ps, dir, frame.size1-1, frame.size2-1);
        const plane = new THREE.Mesh(planeGeometry, material);

        const frameGroup = new THREE.Group();
        frameGroup.add(particles);
        frameGroup.add(plane);
        frames.push(frameGroup);
        scene.add(frameGroup);
        t++;
    }
    currentFrame = 0;
    const velocity = 0.000001;
    const updateFrame = (steps=20) => {
        frames.forEach((f,i) => {
            const dt = processedData[currentFrame].time - processedData[i].time;
            const newPos = dir.clone().multiplyScalar(dt * velocity);
            f.position.lerp(newPos, Math.sqrt(1/steps));
            f.visible = currentFrame >= i;
        });
        render();
        if (steps > 1) {
            requestAnimationFrame(()=>{
                updateFrame(steps-1);
            })
        }
    };
    updateFrame();

    window.addEventListener("keydown", (event) => {
        switch (event.code) {
            case "ArrowRight":
                currentFrame = Math.min(currentFrame+1, frames.length-1);
                updateFrame();
                break;
            case "ArrowLeft":
                currentFrame = Math.max(currentFrame-1, 0);
                updateFrame();
                break;
        }
    });


}

function generateTexture(data, height, width) {
    const canvas = document.createElement("canvas");
    canvas.width = width;
    canvas.height = height;

    const context = canvas.getContext("2d");
    context.fillStyle = "#000";
    context.fillRect(0, 0, width, height);

    const image = context.getImageData(0, 0, canvas.width, canvas.height);
    const imageData = image.data;

    const lut = new Lut("ylOrRd", 512);
    lut.minV = Math.min(...data);
    lut.maxV = Math.max(...data);

    for (let i = 0, j = 0, l = imageData.length; i < l; i += 4, j++) {
        const color = lut.getColor(data[j]);
        imageData[i] = color.r * 255;       // R
        imageData[i + 1] = color.g * 255;   // G
        imageData[i + 2] = color.b * 255;   // B
        imageData[i + 3] = (data[j]/lut.maxV) * 255 * 0.75;   // A
    }

    context.putImageData(image, 0, 0);

    return canvas;
}

function save(blob, filename) {
    const link = document.createElement( 'a' );
    link.style.display = 'none';
    document.body.appendChild(link);

    link.href = URL.createObjectURL(blob);
    link.download = filename;
    link.click();
}

function saveArrayBuffer(buffer, filename) {
    save(new Blob([buffer], {
        type: 'application/octet-stream'
    }), filename);
}

function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
    render();
}

function render() {
    renderer.render(scene, camera);
}