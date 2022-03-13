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


package com.garmin.fit.plugins;

public class ActivityFileValidationResult {

    private final String name;
    private final Level level;
    private String description;
    private Status status;

    public ActivityFileValidationResult(String name, Level level) {
        this.name = name;
        this.level = level;
        this.status = Status.SKIPPED;
    }

    public String getName() {
        return name;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public Status getStatus() {
        return status;
    }

    public void setStatus(Status status) {
        this.status = status;
    }

    public Level getLevel() {
        return level;
    }

    @Override
    public String toString() {
        return getName() + " - Level: " + getLevel() + " Status: " + getStatus();
    }

    public enum Status {
        SKIPPED, WARNING, FAILED, PASSED
    }

    public enum Level {
        REQUIRED, OPTIONAL
    }
}
