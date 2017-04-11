from astrocook import spec1d
#from astropy import units as u
import numpy as np
#import matplotlib.pyplot as plt
import copy




class spec1dContinuum():
    """Class for to handle continuum in spectra
    """

    def findStretchable(self,
                        spec,
                        stiff=None,
                        pnratio=None,
                        minStep=None,
                        tol=None):
        """Find absorbed continuum.

        """

        if (stiff is None):
            stiff = 0.9

        if (pnratio is None):
            pnratio = 20.

        if (minStep is None):
            minStep = 5e-4

        if (tol is None):
            tol = 0.01

        if (abs(stiff) > 1):
            raise ValueError('Stiff must be a number between 0 and 1')

        if (abs(stiff) <= 0):
            raise ValueError('pnRatio must be a number > 0')

        cStr = stiff
        cUp  = 1.-stiff
        cDn  = cUp / float(pnratio)
        print("Settings: ", cStr, cUp, cDn)

        offsetX = spec.x.min().value
        dataX   = spec.x.value - offsetX
        scaleX  = dataX.max()
        dataX  /= scaleX

        offsetY = spec.y.min().value
        dataY   = spec.y.value - offsetY
        scaleY  = dataY.max()
        dataY  /= scaleY

        deltaXsq = ( dataX[1:] - dataX[0:-1] )**2.
        deltaXsq = np.insert(deltaXsq, 0, deltaXsq[0])

        around = 1
        npar = 2
        modelY = None
        lastTotScore = None
        while True:
            npar += (npar-1)

            if (npar > len(dataX)):
                print('Algorithm did not converged')
                break

            index = np.int_(np.round(np.linspace(0, len(dataX)-1, npar)))
            if (modelY is not None):
                newModelX = dataX[index]
                modelY = np.interp(newModelX, modelX, modelY)
                modelX = copy.deepcopy(newModelX)
            else:
                modelX = dataX[index]
                modelY = np.interp(modelX, dataX, dataY)

            go = np.ones(npar)
            cur = -1
            while go.sum() > 0:
                cur = (cur+1) % npar
                if (go[cur] == 0):
                    cur += 1
                    continue

                i0 = (cur-around)
                if (i0 < 0):
                    i0 = 0

                i1 = (cur+around)
                if (i1>(npar-1)):
                    i1 = npar-1

                jj = (dataX >= modelX[i0])  &  (dataX <= modelX[i1])
                valAtStart = modelY[cur]
                step = 0.1
                while step > minStep:
                    for sign in [1, -1]:

                        lastScore = -1
                        while True:
                            if (lastScore >= 0):
                                modelY[cur] += sign * step

                            modelCmp = np.interp(dataX[jj], modelX[i0:i1+1], modelY[i0:i1+1])
                            deltaModelSq = (modelCmp[1:] - modelCmp[0:-1])**2.
                            deltaModelSq = np.insert(deltaModelSq, 0, 0)

                            tmp = dataY[jj]-modelCmp
                            iUp = (tmp > 0)
                            iDn = (tmp < 0)
                            score = \
                                    cStr * np.sqrt(deltaModelSq + deltaXsq[jj]).sum() + \
                                    cUp  * tmp[iUp].sum()                             - \
                                    cDn  * tmp[iDn].sum()

                            if (lastScore != -1):
                                if (score >= lastScore):
                                    modelY[cur] -= sign * step #revert
                                    break

                            lastScore = score
                    step /= 20.

                diff = abs(modelY[cur] - valAtStart)
                if (diff < minStep):
                    go[cur] = 0
                else:
                    go[i0:i1+1] = 1

            modelCmp = np.interp(dataX, modelX, modelY)
            deltaModelSq = (modelCmp[1:] - modelCmp[0:-1])**2.
            deltaModelSq = np.insert(deltaModelSq, 0, 0)
            tmp = dataY-modelCmp
            iUp = (tmp > 0)
            iDn = (tmp < 0)
            score1 = np.sqrt(deltaModelSq + deltaXsq).sum()
            score2 = tmp[iUp].sum()
            score3 = tmp[iDn].sum()
            totScore = cStr * score1 + \
                       cUp  * score2 - \
                       cDn  * score3
            print('Scores: ', score1, score2, score3, totScore, len(modelX))
            if (lastTotScore is not None):
                if (abs(lastTotScore - totScore)/totScore < tol):
                    break

            lastTotScore = totScore

            #plt.plot(dataX, dataY, marker='.', linestyle='None', markersize=1)
            #plt.plot(dataX, modelCmp)
            #plt.show()

        i = np.argwhere(modelCmp > dataY)
        print('Length:   ', score1)
        print('PN ratio: ', len(i)/float(len(dataX) - len(i)), ' ~ ', pnratio)

        modelCmp = np.interp(dataX, modelX, modelY)
        modelCmp *= scaleY
        modelCmp += offsetY

        return modelCmp
