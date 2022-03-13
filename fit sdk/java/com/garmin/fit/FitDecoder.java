////////////////////////////////////////////////////////////////////////////////
// The following FIT Protocol software provided may be used with FIT protocol
// devices only and remains the copyrighted property of Garmin Canada Inc.
// The software is being provided on an "as-is" basis and as an accommodation,
// and therefore all warranties, representations, or guarantees of any kind
// (whether express, implied or statutory) including, without limitation,
// warranties of merchantability, non-infringement, or fitness for a particular
// purpose, are specifically disclaimed.
//
// Copyright 2021 Garmin International, Inc.
////////////////////////////////////////////////////////////////////////////////
// ****WARNING****  This file is auto-generated!  Do NOT edit this file.
// Profile Version = 21.67Release
// Tag = production/akw/21.67.00-0-gd790f76b
////////////////////////////////////////////////////////////////////////////////

package com.garmin.fit;

import com.garmin.fit.BufferedMesgBroadcaster;
import com.garmin.fit.Decode;
import com.garmin.fit.FitRuntimeException;
import com.garmin.fit.MesgBroadcastPlugin;

import java.io.InputStream;

public class FitDecoder {

    public FitDecoder() {
    }

    public FitMessages decode(InputStream inputStream) {
        Decode decode = new Decode();
        FitListener fitListener = new FitListener();

        try {
            decode.read(inputStream, fitListener, null);
        } catch (FitRuntimeException fre) {
            // If a FIT file with 0 data size is encountered, attempt to
            // process the next chained FIT file in the input stream.
            if (decode.getInvalidFileDataSize()) {
                decode.nextFile();
                decode.read(inputStream, fitListener, null);
            } else {
                throw fre;
            }
        }

        return fitListener.getFitMessages();
    }

    public FitMessages decode(InputStream inputStream, MesgBroadcastPlugin plugin) {
        Decode decode = new Decode();
        FitListener fitListener = new FitListener();
        BufferedMesgBroadcaster mesgBroadcaster = new BufferedMesgBroadcaster(decode);

        mesgBroadcaster.registerMesgBroadcastPlugin(plugin);
        mesgBroadcaster.addListener(fitListener);

        try {
            mesgBroadcaster.run(inputStream);
            mesgBroadcaster.broadcast();
        } catch (FitRuntimeException fre) {
            // If a FIT file with 0 data size is encountered, attempt to
            // process the next chained FIT file in the input stream.
            if (decode.getInvalidFileDataSize()) {
                decode.nextFile();
                mesgBroadcaster.run(inputStream);
                mesgBroadcaster.broadcast();
            } else {
                throw fre;
            }
        }

        return fitListener.getFitMessages();
    }
}