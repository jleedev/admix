class CreateGenotypeFiles < ActiveRecord::Migration
  def self.up
    create_table :genotype_files do |t|
      t.belongs_to :user

      t.string :name
      t.text :data

      t.integer :num_alleles
      t.integer :num_rows

      t.timestamps
    end
  end

  def self.down
    drop_table :genotype_files
  end
end
