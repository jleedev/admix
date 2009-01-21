class CreateLocusFiles < ActiveRecord::Migration
  def self.up
    create_table :locus_files do |t|
      t.integer :population
      t.belongs_to :user

      t.timestamps
    end
  end

  def self.down
    drop_table :locus_files
  end
end
