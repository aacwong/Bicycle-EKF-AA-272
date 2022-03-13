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


public enum WeatherStatus  {
    CLEAR((short)0),
    PARTLY_CLOUDY((short)1),
    MOSTLY_CLOUDY((short)2),
    RAIN((short)3),
    SNOW((short)4),
    WINDY((short)5),
    THUNDERSTORMS((short)6),
    WINTRY_MIX((short)7),
    FOG((short)8),
    HAZY((short)11),
    HAIL((short)12),
    SCATTERED_SHOWERS((short)13),
    SCATTERED_THUNDERSTORMS((short)14),
    UNKNOWN_PRECIPITATION((short)15),
    LIGHT_RAIN((short)16),
    HEAVY_RAIN((short)17),
    LIGHT_SNOW((short)18),
    HEAVY_SNOW((short)19),
    LIGHT_RAIN_SNOW((short)20),
    HEAVY_RAIN_SNOW((short)21),
    CLOUDY((short)22),
    INVALID((short)255);

    protected short value;

    private WeatherStatus(short value) {
        this.value = value;
    }

    public static WeatherStatus getByValue(final Short value) {
        for (final WeatherStatus type : WeatherStatus.values()) {
            if (value == type.value)
                return type;
        }

        return WeatherStatus.INVALID;
    }

    /**
     * Retrieves the String Representation of the Value
     * @return The string representation of the value
     */
    public static String getStringFromValue( WeatherStatus value ) {
        return value.name();
    }

    public short getValue() {
        return value;
    }


}
