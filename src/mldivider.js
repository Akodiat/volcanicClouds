/**
 * Divide columns by matrix (matlab mldivide)
 * The code was compiled from matlab into C++ and then translated to
 * JavaScript, so the magic numbers, etc, can indeed be regarded as
 * magic...
 * @param {*} cols
 * @param {*} matrix
 * @returns
 */
function mldivider(cols, matrix) {
    // Flatten matrices
    const [size1, size2] = matrix.size();
    const flatCols = cols.toArray().flatMap(v=>v);
    const flatMatrix = matrix.toArray().flatMap(v=>v);
    let concentration = new Array(size2);

    let tau = 0.0;
    let wj;
    let xnorm_tmp_tmp_tmp;
    let knt;
    let rankA = 0;

    // Perform division
    wj = flatCols[0];
    tau = 0.0;
    xnorm_tmp_tmp_tmp = xnrm2(flatCols);
    if (xnorm_tmp_tmp_tmp !== 0.0) {
        let beta1 = rt_hypotd_snf(flatCols[0], xnorm_tmp_tmp_tmp);
        if (flatCols[0] >= 0.0) {
            beta1 = -beta1;
        }
        if (Math.abs(beta1) < 1.0020841800044864e-292) {
            knt = 0;
            do {
                knt++;
                for (let b_k = 0; b_k < size1-1; b_k++) {
                    flatCols[b_k + 1] *= 9.9792015476736e+291;
                }
                beta1 *= 9.9792015476736e+291;
                wj *= 9.9792015476736e+291;
            } while (Math.abs(beta1) < 1.0020841800044864e-292 && knt < 20);
            beta1 = rt_hypotd_snf(wj, xnrm2(flatCols));
            if (wj >= 0.0) {
                beta1 = -beta1;
            }
            tau = (beta1 - wj) / beta1;
            xnorm_tmp_tmp_tmp = 1.0 / (wj - beta1);
            for (let b_k = 0; b_k < size1-1; b_k++) {
                flatCols[b_k + 1] *= xnorm_tmp_tmp_tmp;
            }
            for (let b_k = 0; b_k < knt; b_k++) {
                beta1 *= 1.0020841800044864e-292;
            }
            wj = beta1;
        } else {
            tau = (beta1 - flatCols[0]) / beta1;
            xnorm_tmp_tmp_tmp = 1.0 / (flatCols[0] - beta1);
            for (let b_k = 0; b_k < size1-1; b_k++) {
                flatCols[b_k + 1] *= xnorm_tmp_tmp_tmp;
            }
            wj = beta1;
        }
    }
    flatCols[0] = wj;

    xnorm_tmp_tmp_tmp = Math.abs(flatCols[0]);
    if (!(xnorm_tmp_tmp_tmp <= 1.4654943925052066e-13 * xnorm_tmp_tmp_tmp)) {
        rankA = 1;
    }

    concentration.fill(0.0);
    if (tau !== 0.0) {
        for (let k = 0; k < size2; k++) {
            xnorm_tmp_tmp_tmp = flatMatrix[size1 * k];
            wj = xnorm_tmp_tmp_tmp;
            for (let b_k = 0; b_k < size1-1; b_k++) {
                wj += flatCols[b_k + 1] * flatMatrix[(b_k + size1 * k)];
            }
            wj *= tau;
            if (wj !== 0.0) {
                flatMatrix[size1 * k] = xnorm_tmp_tmp_tmp - wj;
                for (let b_k = 0; b_k < size1-1; b_k++) {
                    knt = (b_k + size1 * k) + 1;
                    flatMatrix[knt] -= flatCols[b_k+1] * wj;
                }
            }
        }
    }

    for (let k = 0; k < size2; k++) {
        if (rankA - 1 >= 0) {
            concentration[k] = flatMatrix[size1 * k];
        }
        for (let b_k = rankA; b_k >= 1; b_k--) {
            concentration[k] /= flatCols[0];
        }
    }

    return concentration;
}

/**
 * Calculates the 2-norm of a vector.
 * @param {Number[]} x
 * @returns {Number}
 */
function xnrm2(x) {
    /*
    let scale = 0.0;
    let ssq = 1.0;
    for (let i = 0; i < x.length; i++) {
        if (x[i] !== 0.0) {
            const absxi = Math.abs(x[i]);
            if (scale < absxi) {
                ssq = 1.0 + ssq * (scale / absxi) ** 2;
                scale = absxi;
            } else {
                ssq += (absxi / scale) ** 2;
            }
        }
    }
    return scale * Math.sqrt(ssq);
    */
    let sum = 0;
    for (let i = 0; i < x.length; i++) {
        sum += x[i] ** 2;
    }
    return Math.sqrt(sum);
}

/**
 *
 * @param {Number} u0
 * @param {Number} u1
 * @returns
 */
function rt_hypotd_snf(u0, u1) {
    let a;
    let b;
    let y;
    a = Math.abs(u0);
    b = Math.abs(u1);
    if (a < b) {
        a /= b;
        y = b * Math.sqrt(a * a + 1.0);
    } else if (a > b) {
        b /= a;
        y = a * Math.sqrt(b * b + 1.0);
    } else if (isNaN(b)) {
        y = NaN;
    } else {
        y = a * 1.4142135623730951;
    }
    return y;
}

export {mldivider};