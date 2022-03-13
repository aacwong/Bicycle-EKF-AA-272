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


public enum BacklightMode  {
    OFF((short)0),
    MANUAL((short)1),
    KEY_AND_MESSAGES((short)2),
    AUTO_BRIGHTNESS((short)3),
    SMART_NOTIFICATIONS((short)4),
    KEY_AND_MESSAGES_NIGHT((short)5),
    KEY_AND_MESSAGES_AND_SMART_NOTIFICATIONS((short)6),
    INVALID((short)255);

    protected short value;

    private BacklightMode(short value) {
        this.value = value;
    }

    public static BacklightMode getByValue(final Short value) {
        for (final BacklightMode type : BacklightMode.values()) {
            if (value == type.value)
                return type;
        }

        return BacklightMode.INVALID;
    }

    /**
     * Retrieves the String Representation of the Value
     * @return The string representation of the value
     */
    public static String getStringFromValue( BacklightMode value ) {
        return value.name();
    }

    public short getValue() {
        return value;
    }


}
