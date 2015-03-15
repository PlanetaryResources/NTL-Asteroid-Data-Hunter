
-- -----------------------------------------------------
-- Table `detection_item`
-- -----------------------------------------------------

CREATE TABLE `detection_item` (
  `id` BIGINT NOT NULL AUTO_INCREMENT,
  `timestamp` TIMESTAMP NOT NULL,
  `neo` TINYINT(1) NOT NULL,
  `submitted` TINYINT(1) NOT NULL,
  `known_by_mpc` TINYINT(1) NOT NULL,
  `significance` DOUBLE NOT NULL,
  `observatory_code` VARCHAR(100) NOT NULL,
  `image_width` DOUBLE NOT NULL,
  `image_height` DOUBLE NOT NULL,
  `image_id` BIGINT NOT NULL,
  PRIMARY KEY (`id`))
ENGINE = InnoDB charset=utf8;


-- -----------------------------------------------------
-- Table `detection_item_frame`
-- -----------------------------------------------------

CREATE TABLE `detection_item_frame` (
  `detection_item_id` BIGINT NOT NULL,
  `frame` TINYINT NOT NULL,
  `right_ascension` DOUBLE NOT NULL,
  `declination` DOUBLE NOT NULL,
  `rough_magnitude` DOUBLE NOT NULL,
  `observation_latitude` DOUBLE NOT NULL,
  `observation_longitude` DOUBLE NOT NULL,
  `x` INT NOT NULL,
  `y` INT NOT NULL,
  `visualization_image` VARCHAR(255) NOT NULL,
  PRIMARY KEY (`detection_item_id`, `frame`),
  CONSTRAINT `fk_detection_item`
    FOREIGN KEY (`detection_item_id`)
    REFERENCES `detection_item` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB charset=utf8;


-- -----------------------------------------------------
-- Table `help_topic`
-- -----------------------------------------------------

CREATE TABLE `help_topic` (
  `id` BIGINT NOT NULL AUTO_INCREMENT,
  `name` VARCHAR(100) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE INDEX `name_UNIQUE` (`name` ASC))
ENGINE = InnoDB charset=utf8;


-- -----------------------------------------------------
-- Table `help_item`
-- -----------------------------------------------------

CREATE TABLE `help_item` (
  `id` BIGINT NOT NULL AUTO_INCREMENT,
  `help_topic_id` BIGINT NOT NULL,
  `title` VARCHAR(255) NOT NULL,
  `content` TEXT NOT NULL,
  PRIMARY KEY (`id`),
  CONSTRAINT `fk_help_topic`
    FOREIGN KEY (`help_topic_id`)
    REFERENCES `help_topic` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB charset=utf8;

CREATE INDEX fk_help_topic_idx ON help_item(help_topic_id ASC);  


-- -----------------------------------------------------
-- Table `detection_session`
-- -----------------------------------------------------

CREATE TABLE `detection_session` (
  `id` BIGINT NOT NULL AUTO_INCREMENT,
  `progress` FLOAT NOT NULL,
  PRIMARY KEY (`id`))
ENGINE = InnoDB charset=utf8;

-- -----------------------------------------------------
-- Table `image_info`
-- -----------------------------------------------------

CREATE TABLE `image_info` (
  `id` BIGINT NOT NULL AUTO_INCREMENT,
  `time_obs` VARCHAR(255) NOT NULL,
  `date_obs` VARCHAR(255) NOT NULL,
  `ra` VARCHAR(255) NOT NULL,
  `dec` VARCHAR(255) NOT NULL,
  `image_create_date` VARCHAR(255) NOT NULL,
  `timestamp` TIMESTAMP NOT NULL,
  `observatory_code` VARCHAR(100) NOT NULL,
  PRIMARY KEY (`id`))
ENGINE = InnoDB charset=utf8;
