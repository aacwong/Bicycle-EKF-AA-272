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



public class StressLevelMesg extends Mesg   {

    
    public static final int StressLevelValueFieldNum = 0;
    
    public static final int StressLevelTimeFieldNum = 1;
    

    protected static final  Mesg stressLevelMesg;
    static {
        // stress_level
        stressLevelMesg = new Mesg("stress_level", MesgNum.STRESS_LEVEL);
        stressLevelMesg.addField(new Field("stress_level_value", StressLevelValueFieldNum, 131, 1, 0, "", false, Profile.Type.SINT16));
        
        stressLevelMesg.addField(new Field("stress_level_time", StressLevelTimeFieldNum, 134, 1, 0, "s", false, Profile.Type.DATE_TIME));
        
    }

    public StressLevelMesg() {
        super(Factory.createMesg(MesgNum.STRESS_LEVEL));
    }

    public StressLevelMesg(final Mesg mesg) {
        super(mesg);
    }


    /**
     * Get stress_level_value field
     *
     * @return stress_level_value
     */
    public Short getStressLevelValue() {
        return getFieldShortValue(0, 0, Fit.SUBFIELD_INDEX_MAIN_FIELD);
    }

    /**
     * Set stress_level_value field
     *
     * @param stressLevelValue
     */
    public void setStressLevelValue(Short stressLevelValue) {
        setFieldValue(0, 0, stressLevelValue, Fit.SUBFIELD_INDEX_MAIN_FIELD);
    }

    /**
     * Get stress_level_time field
     * Units: s
     * Comment: Time stress score was calculated
     *
     * @return stress_level_time
     */
    public DateTime getStressLevelTime() {
        return timestampToDateTime(getFieldLongValue(1, 0, Fit.SUBFIELD_INDEX_MAIN_FIELD));
    }

    /**
     * Set stress_level_time field
     * Units: s
     * Comment: Time stress score was calculated
     *
     * @param stressLevelTime
     */
    public void setStressLevelTime(DateTime stressLevelTime) {
        setFieldValue(1, 0, stressLevelTime.getTimestamp(), Fit.SUBFIELD_INDEX_MAIN_FIELD);
    }

}
