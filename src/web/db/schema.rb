# This file is auto-generated from the current state of the database. Instead of editing this file, 
# please use the migrations feature of Active Record to incrementally modify your database, and
# then regenerate this schema definition.
#
# Note that this schema.rb definition is the authoritative source for your database schema. If you need
# to create the application database on another system, you should be using db:schema:load, not running
# all the migrations from scratch. The latter is a flawed and unsustainable approach (the more migrations
# you'll amass, the slower it'll run and the greater likelihood for issues).
#
# It's strongly recommended to check this file into your version control system.

ActiveRecord::Schema.define(:version => 20090224065233) do

  create_table "allele_frequencies", :force => true do |t|
    t.float    "freq"
    t.integer  "population_number"
    t.integer  "marker_id"
    t.datetime "created_at"
    t.datetime "updated_at"
  end

  create_table "genotype_files", :force => true do |t|
    t.integer  "user_id"
    t.string   "name"
    t.text     "data"
    t.integer  "num_markers"
    t.integer  "num_rows"
    t.datetime "created_at"
    t.datetime "updated_at"
  end

  create_table "jobs", :force => true do |t|
    t.datetime "created_at"
    t.datetime "updated_at"
  end

  create_table "locus_files", :force => true do |t|
    t.integer  "population"
    t.integer  "user_id"
    t.datetime "created_at"
    t.datetime "updated_at"
    t.string   "name"
  end

  create_table "markers", :force => true do |t|
    t.string   "name",          :limit => 16
    t.integer  "locus_file_id"
    t.datetime "created_at"
    t.datetime "updated_at"
  end

  create_table "users", :force => true do |t|
    t.string   "name"
    t.datetime "created_at"
    t.datetime "updated_at"
  end

end
