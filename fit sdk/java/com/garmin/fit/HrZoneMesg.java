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



public class HrZoneMesg extends Mesg   {

    
    public static final int MessageIndexFieldNum = 254;
    
    public static final int HighBpmFieldNum = 1;
    
    public static final int NameFieldNum = 2;
    

    protected static final  Mesg hrZoneMesg;
    static {
        // hr_zone
        hrZoneMesg = new Mesg("hr_zone", MesgNum.HR_ZONE);
        hrZoneMesg.addField(new Field("message_index", MessageIndexFieldNum, 132, 1, 0, "", false, Profile.Type.MESSAGE_INDEX));
        
        hrZoneMesg.addField(new Field("high_bpm", HighBpmFieldNum, 2, 1, 0, "bpm", false, Profile.Type.UINT8));
        
        hrZoneMesg.addField(new Field("name", NameFieldNum, 7, 1, 0, "", false, Profile.Type.STRING));
        
    }

    public HrZoneMesg() {
        super(Factory.createMesg(MesgNum.HR_ZONE));
    }

    public HrZoneMesg(final Mesg mesg) {
        super(mesg);
    }


    /**
     * Get message_index field
     *
     * @return message_index
     */
    public Integer getMessageIndex() {
        return getFieldIntegerValue(254, 0, Fit.SUBFIELD_INDEX_MAIN_FIELD);
    }

    /**
     * Set message_index field
     *
     * @param messageIndex
     */
    public void setMessageIndex(Integer messageIndex) {
        setFieldValue(254, 0, messageIndex, Fit.SUBFIELD_INDEX_MAIN_FIELD);
    }

    /**
     * Get high_bpm field
     * Units: bpm
     *
     * @return high_bpm
     */
    public Short getHighBpm() {
        return getFieldShortValue(1, 0, Fit.SUBFIELD_INDEX_MAIN_FIELD);
    }

    /**
     * Set high_bpm field
     * Units: bpm
     *
     * @param highBpm
     */
    public void setHighBpm(Short highBpm) {
        setFieldValue(1, 0, highBpm, Fit.SUBFIELD_INDEX_MAIN_FIELD);
    }

    /**
     * Get name field
     *
     * @return name
     */
    public String getName() {
        return getFieldStringValue(2, 0, Fit.SUBFIELD_INDEX_MAIN_FIELD);
    }

    /**
     * Set name field
     *
     * @param name
     */
    public void setName(String name) {
        setFieldValue(2, 0, name, Fit.SUBFIELD_INDEX_MAIN_FIELD);
    }

}
