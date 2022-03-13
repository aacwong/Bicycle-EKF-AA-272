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

import java.util.HashMap;
import java.util.Map;

public class TricepsExtensionExerciseName  {
    public static final int BENCH_DIP = 0;
    public static final int WEIGHTED_BENCH_DIP = 1;
    public static final int BODY_WEIGHT_DIP = 2;
    public static final int CABLE_KICKBACK = 3;
    public static final int CABLE_LYING_TRICEPS_EXTENSION = 4;
    public static final int CABLE_OVERHEAD_TRICEPS_EXTENSION = 5;
    public static final int DUMBBELL_KICKBACK = 6;
    public static final int DUMBBELL_LYING_TRICEPS_EXTENSION = 7;
    public static final int EZ_BAR_OVERHEAD_TRICEPS_EXTENSION = 8;
    public static final int INCLINE_DIP = 9;
    public static final int WEIGHTED_INCLINE_DIP = 10;
    public static final int INCLINE_EZ_BAR_LYING_TRICEPS_EXTENSION = 11;
    public static final int LYING_DUMBBELL_PULLOVER_TO_EXTENSION = 12;
    public static final int LYING_EZ_BAR_TRICEPS_EXTENSION = 13;
    public static final int LYING_TRICEPS_EXTENSION_TO_CLOSE_GRIP_BENCH_PRESS = 14;
    public static final int OVERHEAD_DUMBBELL_TRICEPS_EXTENSION = 15;
    public static final int RECLINING_TRICEPS_PRESS = 16;
    public static final int REVERSE_GRIP_PRESSDOWN = 17;
    public static final int REVERSE_GRIP_TRICEPS_PRESSDOWN = 18;
    public static final int ROPE_PRESSDOWN = 19;
    public static final int SEATED_BARBELL_OVERHEAD_TRICEPS_EXTENSION = 20;
    public static final int SEATED_DUMBBELL_OVERHEAD_TRICEPS_EXTENSION = 21;
    public static final int SEATED_EZ_BAR_OVERHEAD_TRICEPS_EXTENSION = 22;
    public static final int SEATED_SINGLE_ARM_OVERHEAD_DUMBBELL_EXTENSION = 23;
    public static final int SINGLE_ARM_DUMBBELL_OVERHEAD_TRICEPS_EXTENSION = 24;
    public static final int SINGLE_DUMBBELL_SEATED_OVERHEAD_TRICEPS_EXTENSION = 25;
    public static final int SINGLE_LEG_BENCH_DIP_AND_KICK = 26;
    public static final int WEIGHTED_SINGLE_LEG_BENCH_DIP_AND_KICK = 27;
    public static final int SINGLE_LEG_DIP = 28;
    public static final int WEIGHTED_SINGLE_LEG_DIP = 29;
    public static final int STATIC_LYING_TRICEPS_EXTENSION = 30;
    public static final int SUSPENDED_DIP = 31;
    public static final int WEIGHTED_SUSPENDED_DIP = 32;
    public static final int SWISS_BALL_DUMBBELL_LYING_TRICEPS_EXTENSION = 33;
    public static final int SWISS_BALL_EZ_BAR_LYING_TRICEPS_EXTENSION = 34;
    public static final int SWISS_BALL_EZ_BAR_OVERHEAD_TRICEPS_EXTENSION = 35;
    public static final int TABLETOP_DIP = 36;
    public static final int WEIGHTED_TABLETOP_DIP = 37;
    public static final int TRICEPS_EXTENSION_ON_FLOOR = 38;
    public static final int TRICEPS_PRESSDOWN = 39;
    public static final int WEIGHTED_DIP = 40;
    public static final int INVALID = Fit.UINT16_INVALID;

    private static final Map<Integer, String> stringMap;

    static {
        stringMap = new HashMap<Integer, String>();
        stringMap.put(BENCH_DIP, "BENCH_DIP");
        stringMap.put(WEIGHTED_BENCH_DIP, "WEIGHTED_BENCH_DIP");
        stringMap.put(BODY_WEIGHT_DIP, "BODY_WEIGHT_DIP");
        stringMap.put(CABLE_KICKBACK, "CABLE_KICKBACK");
        stringMap.put(CABLE_LYING_TRICEPS_EXTENSION, "CABLE_LYING_TRICEPS_EXTENSION");
        stringMap.put(CABLE_OVERHEAD_TRICEPS_EXTENSION, "CABLE_OVERHEAD_TRICEPS_EXTENSION");
        stringMap.put(DUMBBELL_KICKBACK, "DUMBBELL_KICKBACK");
        stringMap.put(DUMBBELL_LYING_TRICEPS_EXTENSION, "DUMBBELL_LYING_TRICEPS_EXTENSION");
        stringMap.put(EZ_BAR_OVERHEAD_TRICEPS_EXTENSION, "EZ_BAR_OVERHEAD_TRICEPS_EXTENSION");
        stringMap.put(INCLINE_DIP, "INCLINE_DIP");
        stringMap.put(WEIGHTED_INCLINE_DIP, "WEIGHTED_INCLINE_DIP");
        stringMap.put(INCLINE_EZ_BAR_LYING_TRICEPS_EXTENSION, "INCLINE_EZ_BAR_LYING_TRICEPS_EXTENSION");
        stringMap.put(LYING_DUMBBELL_PULLOVER_TO_EXTENSION, "LYING_DUMBBELL_PULLOVER_TO_EXTENSION");
        stringMap.put(LYING_EZ_BAR_TRICEPS_EXTENSION, "LYING_EZ_BAR_TRICEPS_EXTENSION");
        stringMap.put(LYING_TRICEPS_EXTENSION_TO_CLOSE_GRIP_BENCH_PRESS, "LYING_TRICEPS_EXTENSION_TO_CLOSE_GRIP_BENCH_PRESS");
        stringMap.put(OVERHEAD_DUMBBELL_TRICEPS_EXTENSION, "OVERHEAD_DUMBBELL_TRICEPS_EXTENSION");
        stringMap.put(RECLINING_TRICEPS_PRESS, "RECLINING_TRICEPS_PRESS");
        stringMap.put(REVERSE_GRIP_PRESSDOWN, "REVERSE_GRIP_PRESSDOWN");
        stringMap.put(REVERSE_GRIP_TRICEPS_PRESSDOWN, "REVERSE_GRIP_TRICEPS_PRESSDOWN");
        stringMap.put(ROPE_PRESSDOWN, "ROPE_PRESSDOWN");
        stringMap.put(SEATED_BARBELL_OVERHEAD_TRICEPS_EXTENSION, "SEATED_BARBELL_OVERHEAD_TRICEPS_EXTENSION");
        stringMap.put(SEATED_DUMBBELL_OVERHEAD_TRICEPS_EXTENSION, "SEATED_DUMBBELL_OVERHEAD_TRICEPS_EXTENSION");
        stringMap.put(SEATED_EZ_BAR_OVERHEAD_TRICEPS_EXTENSION, "SEATED_EZ_BAR_OVERHEAD_TRICEPS_EXTENSION");
        stringMap.put(SEATED_SINGLE_ARM_OVERHEAD_DUMBBELL_EXTENSION, "SEATED_SINGLE_ARM_OVERHEAD_DUMBBELL_EXTENSION");
        stringMap.put(SINGLE_ARM_DUMBBELL_OVERHEAD_TRICEPS_EXTENSION, "SINGLE_ARM_DUMBBELL_OVERHEAD_TRICEPS_EXTENSION");
        stringMap.put(SINGLE_DUMBBELL_SEATED_OVERHEAD_TRICEPS_EXTENSION, "SINGLE_DUMBBELL_SEATED_OVERHEAD_TRICEPS_EXTENSION");
        stringMap.put(SINGLE_LEG_BENCH_DIP_AND_KICK, "SINGLE_LEG_BENCH_DIP_AND_KICK");
        stringMap.put(WEIGHTED_SINGLE_LEG_BENCH_DIP_AND_KICK, "WEIGHTED_SINGLE_LEG_BENCH_DIP_AND_KICK");
        stringMap.put(SINGLE_LEG_DIP, "SINGLE_LEG_DIP");
        stringMap.put(WEIGHTED_SINGLE_LEG_DIP, "WEIGHTED_SINGLE_LEG_DIP");
        stringMap.put(STATIC_LYING_TRICEPS_EXTENSION, "STATIC_LYING_TRICEPS_EXTENSION");
        stringMap.put(SUSPENDED_DIP, "SUSPENDED_DIP");
        stringMap.put(WEIGHTED_SUSPENDED_DIP, "WEIGHTED_SUSPENDED_DIP");
        stringMap.put(SWISS_BALL_DUMBBELL_LYING_TRICEPS_EXTENSION, "SWISS_BALL_DUMBBELL_LYING_TRICEPS_EXTENSION");
        stringMap.put(SWISS_BALL_EZ_BAR_LYING_TRICEPS_EXTENSION, "SWISS_BALL_EZ_BAR_LYING_TRICEPS_EXTENSION");
        stringMap.put(SWISS_BALL_EZ_BAR_OVERHEAD_TRICEPS_EXTENSION, "SWISS_BALL_EZ_BAR_OVERHEAD_TRICEPS_EXTENSION");
        stringMap.put(TABLETOP_DIP, "TABLETOP_DIP");
        stringMap.put(WEIGHTED_TABLETOP_DIP, "WEIGHTED_TABLETOP_DIP");
        stringMap.put(TRICEPS_EXTENSION_ON_FLOOR, "TRICEPS_EXTENSION_ON_FLOOR");
        stringMap.put(TRICEPS_PRESSDOWN, "TRICEPS_PRESSDOWN");
        stringMap.put(WEIGHTED_DIP, "WEIGHTED_DIP");
    }


    /**
     * Retrieves the String Representation of the Value
     * @return The string representation of the value, or empty if unknown
     */
    public static String getStringFromValue( Integer value ) {
        if( stringMap.containsKey( value ) ) {
            return stringMap.get( value );
        }

        return "";
    }

    /**
     * Retrieves a value given a string representation
     * @return The value or INVALID if unkwown
     */
    public static Integer getValueFromString( String value ) {
        for( Map.Entry<Integer, String> entry : stringMap.entrySet() ) {
            if( entry.getValue().equals( value ) ) {
                return entry.getKey();
            }
        }

        return INVALID;
    }

}
